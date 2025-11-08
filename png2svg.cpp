/*png2svg.cpp
// Ultra-fast tool to convert PNG images to SVG
// True raster-to-vector conversion: creates SVG with path data, not embedded bitmaps. 
// Requires libpng: sudo apt update && sudo apt install libpng-dev
// Optional pkg-config: sudo apt update && sudo apt install pkg-config
// Compile: g++ -O3 -std=c++17 -pthread -lpng png2svg.cpp -o png2svg
// Compile with pkg-config: g++ png2svg.cpp -o png2svg $(pkg-config --cflags --libs libpng)
// Tested on Debian linux with libpng
*/

#include <png.h>
#include <bits/stdc++.h>
using namespace std;

// ---------------------- Utilities ----------------------
// Pixel integer point
struct PixelPt { int x,y; };
inline bool operator==(PixelPt const& A, PixelPt const& B){ return A.x==B.x && A.y==B.y; }
inline bool operator!=(PixelPt const& A, PixelPt const& B){ return !(A==B); }

// Pixel integer edge
struct PixelEdge { PixelPt a,b; };
inline bool operator==(PixelEdge const& A, PixelEdge const& B){ return A.a==B.a && A.b==B.b; }

// hash for PixelEdge
struct PixelEdgeHash {
    size_t operator()(PixelEdge const& e) const noexcept {
        // combine coordinates into size_t
        size_t h = 1469598103934665603ULL;
        auto mix = [&](long long v){
            h ^= (size_t)v;
            h *= 1099511628211ULL;
        };
        mix(e.a.x); mix(e.a.y); mix(e.b.x); mix(e.b.y);
        return h;
    }
};

struct Options {
    bool contiguous = true; // default: join contiguous areas
    bool opaque = false;
    bool keep_every_point = false;
    int tile_size = 512;
    int threads = thread::hardware_concurrency() ? (int)thread::hardware_concurrency() : 4;
    string filename;
};

// ---------------------- PNG loader (libpng) ----------------------
bool load_png_rgba(const string &filename, int &width, int &height, vector<unsigned char> &out_rgba){
    FILE *fp = fopen(filename.c_str(), "rb");
    if(!fp) return false;
    png_structp png_ptr = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if(!png_ptr){ fclose(fp); return false; }
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if(!info_ptr){ png_destroy_read_struct(&png_ptr,(png_infopp)0,(png_infopp)0); fclose(fp); return false; }
    if(setjmp(png_jmpbuf(png_ptr))){ png_destroy_read_struct(&png_ptr,&info_ptr,(png_infopp)0); fclose(fp); return false; }
    png_init_io(png_ptr, fp);
    png_read_info(png_ptr, info_ptr);
    width = png_get_image_width(png_ptr, info_ptr);
    height = png_get_image_height(png_ptr, info_ptr);
    png_byte color_type = png_get_color_type(png_ptr, info_ptr);
    png_byte bit_depth  = png_get_bit_depth(png_ptr, info_ptr);
    // conversions to 8-bit RGBA
    if(bit_depth == 16) png_set_strip_16(png_ptr);
    if(color_type == PNG_COLOR_TYPE_PALETTE) png_set_palette_to_rgb(png_ptr);
    if(color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8) png_set_expand_gray_1_2_4_to_8(png_ptr);
    if(png_get_valid(png_ptr, info_ptr, PNG_INFO_tRNS)) png_set_tRNS_to_alpha(png_ptr);
    if(color_type == PNG_COLOR_TYPE_RGB || color_type == PNG_COLOR_TYPE_GRAY) png_set_add_alpha(png_ptr, 0xFF, PNG_FILLER_AFTER);
    if(color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_GRAY_ALPHA) png_set_gray_to_rgb(png_ptr);
    png_read_update_info(png_ptr, info_ptr);

    vector<png_bytep> rows(height);
    vector<unsigned char> img;
    img.resize((size_t)width * height * 4);
    for(int y=0;y<height;++y) rows[y] = (png_bytep)(img.data() + (size_t)y * width * 4);
    png_read_image(png_ptr, rows.data());
    png_destroy_read_struct(&png_ptr, &info_ptr, nullptr);
    fclose(fp);
    out_rgba.swap(img);
    return true;
}

// pack/unpack convenience
inline uint32_t packRGBA(const unsigned char* p){
    return ((uint32_t)p[0]<<24) | ((uint32_t)p[1]<<16) | ((uint32_t)p[2]<<8) | ((uint32_t)p[3]);
}
inline void unpackRGBA(uint32_t v, int &r,int &g,int &b,int &a){
    r = (v>>24)&0xFF; g=(v>>16)&0xFF; b=(v>>8)&0xFF; a=v&0xFF;
}

// ---------------------- Disjoint-set (Union-Find) ----------------------
struct DisjointSet {
    vector<int> parent;
    DisjointSet(int n=0){ parent.reserve(n); for(int i=0;i<n;++i) parent.push_back(i); }
    int add(){ int id = parent.size(); parent.push_back(id); return id; }
    int find(int x){ while(parent[x]!=x){ parent[x]=parent[parent[x]]; x=parent[x]; } return x; }
    void unify(int a,int b){ a=find(a); b=find(b); if(a!=b) parent[b]=a; }
    int size() const { return (int)parent.size(); }
};

// ---------------------- Per-tile labeling ----------------------
struct TileResult {
    int tile_x, tile_y; // tile indices
    int w,h; // tile width/height (may be smaller at edges)
    vector<int> local_labels; // w*h ints
    vector<vector<PixelPt>> label_pixels; // index 1..L, store global coords
    vector<pair<PixelPt,int>> top_border, right_border, bottom_border, left_border;
};

TileResult process_tile_label(const vector<uint32_t> &img, int img_w, int img_h,
                              int tile_x, int tile_y, int tile_size,
                              const Options &opt)
{
    int tx = tile_x * tile_size;
    int ty = tile_y * tile_size;
    int tw = min(tile_size, img_w - tx);
    int th = min(tile_size, img_h - ty);
    TileResult res;
    res.tile_x = tile_x; res.tile_y = tile_y; res.w = tw; res.h = th;
    res.local_labels.assign(tw*th, 0);

    int next_label = 1;
    unordered_map<int,int> eqmap; // local union map
    auto find_local = [&](int a)->int{
        int x=a; while(eqmap.count(x) && eqmap[x]!=x) x=eqmap[x];
        return x;
    };
    auto unify_local = [&](int a,int b){
        int fa=find_local(a), fb=find_local(b);
        if(fa==0 || fb==0) return;
        if(fa!=fb) eqmap[fb]=fa;
    };

    auto get_global_idx = [&](int lx,int ly)->int{ return (ty+ly)*img_w + (tx+lx); };

    for(int y=0;y<th;++y){
        for(int x=0;x<tw;++x){
            int gidx = get_global_idx(x,y);
            uint32_t rgba = img[gidx];
            if(opt.opaque && ((rgba & 0xFF) == 0)) continue;
            int label = 0;
            if(x>0){
                int left_label = res.local_labels[y*tw + (x-1)];
                if(left_label){
                    int gidx_left = get_global_idx(x-1,y);
                    if(img[gidx_left] == rgba) label = left_label;
                }
            }
            if(y>0){
                int top_label = res.local_labels[(y-1)*tw + x];
                if(top_label){
                    int gidx_top = get_global_idx(x,y-1);
                    if(img[gidx_top] == rgba){
                        if(label==0) label = top_label;
                        else if(label != top_label) {
                            if(!eqmap.count(label)) eqmap[label]=label;
                            if(!eqmap.count(top_label)) eqmap[top_label]=top_label;
                            unify_local(label, top_label);
                        }
                    }
                }
            }
            if(label==0){
                label = next_label++;
                eqmap[label]=label;
            }
            res.local_labels[y*tw + x] = label;
        }
    }

    unordered_map<int,int> canonical;
    int canon_count = 0;
    res.label_pixels.resize(next_label+1);
    for(int y=0;y<th;++y){
        for(int x=0;x<tw;++x){
            int l = res.local_labels[y*tw + x];
            if(l==0) continue;
            int root = find_local(l);
            if(!canonical.count(root)) canonical[root] = ++canon_count;
            int glabel = canonical[root];
            res.local_labels[y*tw + x] = glabel;
            res.label_pixels[glabel].push_back(PixelPt{tx + x, ty + y});
        }
    }
    res.label_pixels.resize(canon_count+1);

    for(int x=0;x<tw;++x){
        int l = res.local_labels[0*tw + x]; if(l) res.top_border.emplace_back(PixelPt{tx + x, ty + 0}, l);
        int lb = res.local_labels[(th-1)*tw + x]; if(lb) res.bottom_border.emplace_back(PixelPt{tx + x, ty + (th-1)}, lb);
    }
    for(int y=0;y<th;++y){
        int l = res.local_labels[y*tw + (tw-1)]; if(l) res.right_border.emplace_back(PixelPt{tx + (tw-1), ty + y}, l);
        int ll = res.local_labels[y*tw + 0]; if(ll) res.left_border.emplace_back(PixelPt{tx + 0, ty + y}, ll);
    }
    return res;
}

// ---------------------- Merge tile borders to global components ----------------------
struct BorderEntry { PixelPt p; int tile_index; int local_label; uint32_t color; };

void merge_tiles_and_build_components(
    const vector<TileResult> &tiles,
    int tiles_x, int tiles_y,
    int img_w, int img_h,
    const vector<uint32_t> &img,
    DisjointSet &ds,
    vector<int> &label_to_component // mapping: tile_index_locallabel_id -> global component id
){
    auto make_uid = [&](int tile_idx,int local_label)->uint64_t{
        return ((uint64_t)tile_idx << 32) | (uint32_t)local_label;
    };
    unordered_map<uint64_t,int> uid_to_ds; uid_to_ds.reserve(100000);

    auto ensure_uid = [&](uint64_t uid)->int{
        auto it = uid_to_ds.find(uid);
        if(it!=uid_to_ds.end()) return it->second;
        int nid = ds.add();
        uid_to_ds[uid] = nid;
        return nid;
    };

    auto color_at = [&](int gx,int gy)->uint32_t{
        return img[gy*img_w + gx];
    };

    int tile_count = (int)tiles.size();
    for(int ti=0; ti<tile_count; ++ti){
        const TileResult &t = tiles[ti];
        int tx = t.tile_x, ty = t.tile_y;
        int tw = t.w, th = t.h;
        if(tx + 1 < tiles_x){
            int neighbor_idx = (ty * tiles_x) + (tx + 1);
            const TileResult &nr = tiles[neighbor_idx];
            unordered_map<long long, pair<int,uint32_t>> neighbor_left_map;
            for(auto &be : nr.left_border){
                long long key = ((long long)be.first.x << 32) | (unsigned int)be.first.y;
                uint32_t col = color_at(be.first.x, be.first.y);
                neighbor_left_map[key] = { be.second, col };
            }
            for(auto &be : t.right_border){
                long long key = ((long long)be.first.x << 32) | (unsigned int)be.first.y;
                auto it = neighbor_left_map.find(key);
                if(it != neighbor_left_map.end()){
                    uint32_t c1 = color_at(be.first.x, be.first.y);
                    uint32_t c2 = it->second.second;
                    if(c1 == c2){
                        uint64_t uid1 = make_uid((int)ti, be.second);
                        uint64_t uid2 = make_uid(neighbor_idx, it->second.first);
                        int n1 = ensure_uid(uid1);
                        int n2 = ensure_uid(uid2);
                        ds.unify(n1,n2);
                    }
                }
            }
        }
        if(ty + 1 < tiles_y){
            int neighbor_idx = ((ty + 1) * tiles_x) + tx;
            const TileResult &nb = tiles[neighbor_idx];
            unordered_map<long long, pair<int,uint32_t>> neighbor_top_map;
            for(auto &be : nb.top_border){
                long long key = ((long long)be.first.x << 32) | (unsigned int)be.first.y;
                neighbor_top_map[key] = { be.second, color_at(be.first.x, be.first.y) };
            }
            for(auto &be : t.bottom_border){
                long long key = ((long long)be.first.x << 32) | (unsigned int)be.first.y;
                auto it = neighbor_top_map.find(key);
                if(it != neighbor_top_map.end()){
                    uint32_t c1 = color_at(be.first.x, be.first.y);
                    uint32_t c2 = it->second.second;
                    if(c1 == c2){
                        uint64_t uid1 = make_uid((int)ti, be.second);
                        uint64_t uid2 = make_uid(neighbor_idx, it->second.first);
                        int n1 = ensure_uid(uid1);
                        int n2 = ensure_uid(uid2);
                        ds.unify(n1,n2);
                    }
                }
            }
        }
    }

    // map tile local labels to components
    size_t sum_size = 0;
    vector<int> tile_base_index(tiles.size());
    for(size_t i=0;i<tiles.size();++i){
        tile_base_index[i] = (int)sum_size;
        sum_size += tiles[i].label_pixels.size();
    }
    label_to_component.assign(sum_size, -1);

    unordered_map<int,int> dsroot_to_component;
    int comp_count = 0;
    for(size_t ti=0; ti<tiles.size(); ++ti){
        const TileResult &t = tiles[ti];
        for(size_t local = 0; local < t.label_pixels.size(); ++local){
            if(t.label_pixels[local].empty()) continue;
            uint64_t uid = ((uint64_t)ti << 32) | (uint32_t)local;
            auto it = uid_to_ds.find(uid);
            int comp_id;
            if(it == uid_to_ds.end()){
                int nid = ds.add();
                uid_to_ds[uid] = nid;
                comp_id = nid;
            } else {
                comp_id = ds.find(it->second);
            }
            int root = ds.find(comp_id);
            auto it2 = dsroot_to_component.find(root);
            if(it2 == dsroot_to_component.end()){
                dsroot_to_component[root] = comp_count++;
            }
            int final_comp = dsroot_to_component[root];
            label_to_component[tile_base_index[ti] + (int)local] = final_comp;
        }
    }
}

// ---------------------- Edge extraction (per component) ----------------------
static const array<pair<PixelPt, pair<PixelPt,PixelPt>>,4> EDGE_MAP = {
    make_pair(PixelPt{-1,0}, make_pair(PixelPt{0,0}, PixelPt{0,1})),
    make_pair(PixelPt{0,1},  make_pair(PixelPt{0,1}, PixelPt{1,1})),
    make_pair(PixelPt{1,0},  make_pair(PixelPt{1,1}, PixelPt{1,0})),
    make_pair(PixelPt{0,-1}, make_pair(PixelPt{1,0}, PixelPt{0,0}))
};

void compute_component_edges(const vector<vector<PixelPt>> &component_pixels,
                             int img_w, int img_h,
                             vector<vector<PixelEdge>> &out_component_edges)
{
    out_component_edges.resize(component_pixels.size());
    for(size_t ci=0; ci<component_pixels.size(); ++ci){
        const auto &pixels = component_pixels[ci];
        if(pixels.empty()) continue;
        unordered_set<long long> pixelset;
        pixelset.reserve(pixels.size()*2);
        for(auto &p : pixels) pixelset.insert( ((long long)p.x<<32) | (unsigned int)p.y );
        unordered_set<long long> edgeset;
        vector<PixelEdge> edges_out;
        edges_out.reserve(pixels.size()*4);
        for(auto &coord : pixels){
            for(auto &em : EDGE_MAP){
                PixelPt offset = em.first;
                int nx = coord.x + offset.x, ny = coord.y + offset.y;
                long long nkey = ((long long)nx<<32) | (unsigned int)ny;
                if(pixelset.find(nkey) != pixelset.end()) continue;
                PixelPt start = { coord.x + em.second.first.x, coord.y + em.second.first.y };
                PixelPt end   = { coord.x + em.second.second.x, coord.y + em.second.second.y };
                long long key = ((long long)start.x<<48) ^ ((long long)start.y<<32) ^ ((long long)end.x<<16) ^ (end.y & 0xFFFF);
                if(edgeset.insert(key).second){
                    edges_out.push_back(PixelEdge{start,end});
                }
            }
        }
        out_component_edges[ci] = move(edges_out);
    }
}

// ---------------------- Join edges (approx. original algorithm) ----------------------
// This function now takes a vector<PixelEdge> and returns vector<vector<PixelEdge>> (pieces)
vector<vector<PixelEdge>> join_edges_set(const vector<PixelEdge> &edges, bool keep_every_point){
    unordered_set<PixelEdge, PixelEdgeHash> assorted(edges.begin(), edges.end());
    vector<vector<PixelEdge>> pieces;
    // Directions order: (0,1),(1,0),(0,-1),(-1,0)
    deque<PixelPt> directions = { PixelPt{0,1}, PixelPt{1,0}, PixelPt{0,-1}, PixelPt{-1,0} };
    auto dir_of = [&](const PixelEdge &e)->PixelPt{ return PixelPt{ e.b.x - e.a.x, e.b.y - e.a.y }; };
    auto normalize_dir = [&](PixelPt d)->PixelPt{
        if(d.x>0) return PixelPt{1,0};
        if(d.x<0) return PixelPt{-1,0};
        if(d.y>0) return PixelPt{0,1};
        if(d.y<0) return PixelPt{0,-1};
        return PixelPt{0,0};
    };

    while(!assorted.empty()){
        PixelEdge start = *assorted.begin();
        assorted.erase(assorted.begin());
        vector<PixelEdge> piece;
        piece.push_back(start);
        while(true){
            PixelPt curr_dir = normalize_dir(dir_of(piece.back()));
            // rotate deque until directions[2] == curr_dir (as in original)
            int rc=0;
            while(!(directions[2]==curr_dir) && rc<4){
                PixelPt back = directions.back(); directions.pop_back(); directions.push_front(back); ++rc;
            }
            bool found=false;
            for(int i=1;i<=3;i++){
                PixelPt nd = directions[i];
                PixelEdge next_edge = PixelEdge{ piece.back().b, PixelPt{ piece.back().b.x + nd.x, piece.back().b.y + nd.y } };
                auto it = assorted.find(next_edge);
                if(it != assorted.end()){
                    PixelEdge ne = *it;
                    assorted.erase(it);
                    if(i==2 && !keep_every_point){
                        piece.back().b = ne.b;
                    } else {
                        piece.push_back(ne);
                    }
                    if(piece.front().a == piece.back().b){
                        // closed
                        PixelPt n0 = normalize_dir(dir_of(piece.front()));
                        PixelPt nl = normalize_dir(dir_of(piece.back()));
                        if(!keep_every_point && n0.x==nl.x && n0.y==nl.y){
                            // merge front and back similar to python
                            PixelPt new_b = piece.front().b;
                            piece.back().b = new_b;
                            piece.erase(piece.begin());
                        }
                        pieces.push_back(piece);
                        piece.clear();
                    }
                    found=true;
                    break;
                }
            }
            if(!found){
                if(!piece.empty()) pieces.push_back(piece);
                break;
            }
            if(piece.empty()) break;
        }
    }
    return pieces;
}

// ---------------------- SVG writing ----------------------
string svg_header(int w,int h){
    char buf[512];
    snprintf(buf, sizeof(buf),
"<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n"
"<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \n"
"  \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
"<svg width=\"%d\" height=\"%d\"\n"
"     xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">\n", w,h);
    return string(buf);
}

// ---------------------- Main flow ----------------------
int main(int argc,char** argv){
    Options opt;
    vector<string> args;
    for(int i=1;i<argc;++i){
        string s = argv[i];
        if(s=="-p" || s=="--pixels") opt.contiguous=false;
        else if(s=="-o" || s=="--opaque") opt.opaque=true;
        else if(s=="-1" || s=="--one") opt.keep_every_point=true;
        else if(s=="-t" || s=="--tile"){ if(i+1<argc){ opt.tile_size = atoi(argv[++i]); } }
        else if(s=="-j" || s=="--jobs"){ if(i+1<argc){ opt.threads = atoi(argv[++i]); } }
        else args.push_back(s);
    }
    if(args.size()!=1){
        cerr << "Usage: " << argv[0] << " [options] FILE.png\n"
             << "  -p, --pixels    produce per-pixel rects (fast)\n"
             << "  -o, --opaque    skip fully transparent pixels\n"
             << "  -1, --one       keep intermediate points on straight edges\n"
             << "  -t TILE_SIZE    tile size (default 512)\n"
             << "  -j THREADS      number of worker threads (default hardware concurrency)\n";
        return 1;
    }
    opt.filename = args[0];

    int w,h;
    vector<unsigned char> raw;
    if(!load_png_rgba(opt.filename, w,h, raw)){
        cerr << "Could not open image: " << opt.filename << "\n";
        return 1;
    }
    vector<uint32_t> img(w*h);
    for(int i=0;i<w*h;i++){
        img[i] = packRGBA(&raw[i*4]);
    }

    if(!opt.contiguous){
        cout << svg_header(w,h);
        int tiles_x = (w + opt.tile_size - 1) / opt.tile_size;
        int tiles_y = (h + opt.tile_size - 1) / opt.tile_size;
        mutex out_m;
        vector<thread> workers;
        atomic<int> next_tile(0);
        auto worker = [&](){
            while(true){
                int idx = next_tile.fetch_add(1);
                if(idx >= tiles_x*tiles_y) break;
                int tx = idx % tiles_x;
                int ty = idx / tiles_x;
                int startx = tx*opt.tile_size;
                int starty = ty*opt.tile_size;
                int tw = min(opt.tile_size, w - startx);
                int th = min(opt.tile_size, h - starty);
                stringstream ss;
                ss.setf(std::ios::fixed); ss<<setprecision(3);
                for(int yy=0; yy<th; ++yy){
                    for(int xx=0; xx<tw; ++xx){
                        int gx = startx + xx, gy = starty + yy;
                        uint32_t c = img[gy*w + gx];
                        unsigned char a = (unsigned char)(c & 0xFF);
                        if(opt.opaque && a==0) continue;
                        int r=(c>>24)&0xFF, g=(c>>16)&0xFF, b=(c>>8)&0xFF;
                        double fa = (double)a/255.0;
                        ss << "  <rect x=\"" << gx << "\" y=\"" << gy << "\" width=\"1\" height=\"1\" style=\"fill:rgb("<<r<<","<<g<<","<<b<<"); fill-opacity:" << fa << "; stroke:none;\" />\n";
                    }
                }
                {
                    lock_guard<mutex> lk(out_m);
                    cout << ss.str();
                }
            }
        };
        for(int i=0;i<opt.threads;++i) workers.emplace_back(worker);
        for(auto &t : workers) t.join();
        cout << "</svg>\n";
        return 0;
    }

    int tiles_x = (w + opt.tile_size - 1) / opt.tile_size;
    int tiles_y = (h + opt.tile_size - 1) / opt.tile_size;
    int tile_count = tiles_x * tiles_y;
    vector<TileResult> tiles(tile_count);

    vector<thread> pool;
    atomic<int> next_tile(0);
    auto labeller = [&](){
        while(true){
            int idx = next_tile.fetch_add(1);
            if(idx >= tile_count) break;
            int tx = idx % tiles_x, ty = idx / tiles_x;
            TileResult tr = process_tile_label(img, w, h, tx, ty, opt.tile_size, opt);
            tiles[idx] = move(tr);
        }
    };
    for(int i=0;i<opt.threads;++i) pool.emplace_back(labeller);
    for(auto &t : pool) t.join();

    DisjointSet ds;
    vector<int> label_to_component;
    merge_tiles_and_build_components(tiles, tiles_x, tiles_y, w, h, img, ds, label_to_component);

    vector<int> tile_base(tile_count);
    int base_sum = 0;
    for(int i=0;i<tile_count;++i){
        tile_base[i] = base_sum;
        base_sum += (int)tiles[i].label_pixels.size();
    }
    int comp_count = 0;
    for(int v : label_to_component) if(v>=0) comp_count = max(comp_count, v+1);
    vector<vector<PixelPt>> components(comp_count);
    for(int ti=0; ti<tile_count; ++ti){
        const TileResult &t = tiles[ti];
        for(int loc=0; loc < (int)t.label_pixels.size(); ++loc){
            int mapidx = tile_base[ti] + loc;
            if(mapidx<0 || mapidx >= (int)label_to_component.size()) continue;
            int cid = label_to_component[mapidx];
            if(cid < 0) continue;
            auto &vec = t.label_pixels[loc];
            auto &dest = components[cid];
            dest.insert(dest.end(), vec.begin(), vec.end());
        }
    }

    // Compute edges per component (parallel)
    vector<vector<PixelEdge>> comp_edges;
    comp_edges.resize(components.size());
    atomic<size_t> next_comp(0);
    vector<thread> edge_pool;
    auto edge_worker = [&](){
        while(true){
            size_t ci = next_comp.fetch_add(1);
            if(ci >= components.size()) break;
            vector<vector<PixelEdge>> one;
            compute_component_edges({components[ci]}, w, h, one);
            if(!one.empty()) comp_edges[ci] = move(one[0]);
        }
    };
    for(int i=0;i<opt.threads;++i) edge_pool.emplace_back(edge_worker);
    for(auto &t : edge_pool) t.join();

    // Join edges per component
    vector<vector<vector<PixelEdge>>> joined_all(components.size());
    for(size_t ci=0; ci<components.size(); ++ci){
        if(comp_edges[ci].empty()) continue;
        joined_all[ci] = join_edges_set(comp_edges[ci], opt.keep_every_point);
    }

    // Output SVG
    cout << svg_header(w,h);
    for(size_t ci=0; ci<components.size(); ++ci){
        if(components[ci].empty()) continue;
        PixelPt sample = components[ci][0];
        uint32_t c = img[sample.y * w + sample.x];
        int r,g,b,a; unpackRGBA(c,r,g,b,a);
        double fa = (double)a/255.0;
        if(joined_all[ci].empty()) continue;
        auto &pieces = joined_all[ci];
        for(auto &piece : pieces){
            if(piece.empty()) continue;
            cout << " <path d=\" ";
            PixelPt st = piece.front().a; cout << "M " << st.x << "," << st.y << " ";
            for(auto &e : piece){
                PixelPt here = e.a;
                cout << "L " << here.x << "," << here.y << " ";
            }
            cout << "Z ";
            cout << "\" style=\"fill:rgb("<<r<<","<<g<<","<<b<<"); fill-opacity:" << fixed << setprecision(3) << fa << "; stroke:none;\" />\n";
        }
    }
    cout << "</svg>\n";
    return 0;
}
