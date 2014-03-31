//boost
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

//stl
#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <time.h>
#include <iomanip>

#ifdef MSVC
#include <windows.h>
#endif

#include "Pair.h"
using namespace boost;

typedef adjacency_list<vecS,vecS,undirectedS,property<vertex_name_t,std::string>,property<edge_weight_t, double> > UGraph;
typedef std::map<std::pair<int,int>,std::vector<int> > Paths;
std::ofstream MSSPF_LOG;
//boost random generator
boost::random::mt19937 rand_double_gen(time(0));

void remove_failed_paths( UGraph & G, std::vector<graph_traits<UGraph>::edge_descriptor> & failed_edges_v )
{
    for (auto i = failed_edges_v.begin(); i != failed_edges_v.end(); ++i)
    {
        remove_edge(source(*i,G),target(*i,G),G);
    }
}

// -1: file error.
// -2: format error.
int import_network_list( const std::string & file_path, UGraph & G)
{
    std::ifstream in(file_path);
    if (!in.is_open())
    {
        return -1;
    }

    std::string word_buffer;
    int node_num = 0, edge_num = 0;
    in >> word_buffer;
    if (word_buffer != "Node_num:")
    {
        return -2;
    }
    in >> word_buffer;
    try
    {
        node_num = boost::lexical_cast<int>(word_buffer);
    }
    catch(std::exception & e)
    {
        std::cerr << e.what() << std::endl;
        return -2;
    }

    in >> word_buffer;
    if (word_buffer != "Edge_num:")
    {
        return -2;
    }
    in >> word_buffer;
    try
    {
        edge_num = boost::lexical_cast<int>(word_buffer);
    }
    catch(std::exception & e)
    {
        std::cerr << e.what() <<std::endl;
        return -2;
    }
    std::getline(in,word_buffer);
    std::getline(in,word_buffer);
    G = UGraph(node_num);

    int link_index = 0, source = 0, destination = 0, length = 0,bandwith = 0, propagation_delay = 0;
    for (int i = 0; i < edge_num; ++i)
    {
        in >> link_index >> source >> destination >> length >> bandwith >> propagation_delay;
        boost::random::uniform_real_distribution<> rand_range_double(1,1.0001);
        double r = rand_range_double(rand_double_gen);
//		std::cout << "#rand\t" << r << std::endl;
        add_edge(source,destination,r,G);
    }
    return 0;
}

// A function to randomly select k items from stream[0..n-1].
template <typename T>
std::vector<T> selectKItems(std::vector<T> & stream, int n, int k)
{
    int i;  // index for elements in stream[]

    // reservoir[] is the output array. Initialize it with
    // first k elements from stream[]
    std::vector<T> reservoir(k);
    for (i = 0; i < k; i++)
        reservoir[i] = stream[i];

    // Use a different seed value so that we don't get
    // same result each time we run this program

    // Iterate from the (k+1)th element to nth element
    for (; i < n; i++)
    {
        // Pick a random index from 0 to i.
        int j = rand() % (i+1);

        // If the randomly  picked index is smaller than k, then replace
        // the element present at the index with new element from stream
        if (j < k)
            reservoir[j] = stream[i];
    }
    return reservoir;
}

void search_affected(
    UGraph & G,
    Paths & original_paths,
    std::vector<graph_traits<UGraph>::edge_descriptor> & failed_edges_v,
    std::map<std::pair<int,int>,std::pair<int,int>> & affected_S_D_paths )
{
    for (auto s_d_path = original_paths.begin(); s_d_path != original_paths.end(); ++s_d_path)
    {
        //traversal every edge in the path
        std::vector<int> & path = s_d_path->second;
        int src = s_d_path->first.first;
        int dst = s_d_path->first.second;
        for (auto v_in_path = path.begin(); v_in_path != path.end() - 1; ++v_in_path)
        {
            //traversal each edge in the failed edges vector to see if they match.
            int src_v = *v_in_path;
            int dst_v = *(v_in_path + 1);
            graph_traits<UGraph>::edge_descriptor each_edge = edge(src_v,dst_v,G).first;

            for (auto each_failed_edge = failed_edges_v.begin(); each_failed_edge != failed_edges_v.end(); ++ each_failed_edge)
            {
                //find a path that need to be reroute.
                if(Z::peer_pair(each_failed_edge->m_source,each_failed_edge->m_target)
                        == Z::peer_pair(each_edge.m_source,each_edge.m_target))
                {

                    affected_S_D_paths[std::make_pair(src,dst)] = std::make_pair(src_v,dst_v);
                }
            }
        }
    }

}

bool is_unique(
    std::vector<std::vector<graph_traits<UGraph>::edge_descriptor>> & failed_edges_v,
    std::vector<graph_traits<UGraph>::edge_descriptor> & combo )
{
    bool r = true;
    for (auto i = failed_edges_v.begin(); i != failed_edges_v.end(); ++i)
    {
        std::vector<graph_traits<UGraph>::edge_descriptor> & a_combo = *i;
        std::set<graph_traits<UGraph>::edge_descriptor> combo_set(combo.begin(),combo.end());
        bool combo_equal = true;
        for (auto j = a_combo.begin(); j != a_combo.end(); ++j)
        {
            if (combo_set.find(*j) == combo_set.end())
            {
                combo_equal = false;
                break;
            }
        }
        if (combo_equal == true)
        {
            r = false;
            break;
        }
    }
    return r;
}

bool connected( std::vector<graph_traits<UGraph>::edge_descriptor> combo, const UGraph & G )
{
    UGraph tested_G = G;
    remove_failed_paths(tested_G,combo);
    std::vector<int> component(num_vertices(tested_G));
    int num_comp = connected_components(tested_G,make_iterator_property_map(component.begin(),get(vertex_index,tested_G)));
    if (num_comp == 2)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void select_failed_paths( int L, int times, UGraph & G,
                          std::vector<std::vector<graph_traits<UGraph>::edge_descriptor>> & failed_edges_v)
{
    srand((unsigned int)time(NULL));
    std::vector<graph_traits<UGraph>::edge_descriptor> edges_G(edges(G).first,edges(G).second);
    int max_do = 1000;
    do
    {
        std::vector<graph_traits<UGraph>::edge_descriptor> combo = selectKItems(edges_G,edges_G.size(),L);
        if ( is_unique(failed_edges_v,combo) && connected(combo,G))
        {
            failed_edges_v.push_back(combo);

        }
        --max_do;
    }
    while (failed_edges_v.size() != times && max_do > 0);
}


void shortest_path_s_d( int s, int d, std::vector<int> &shortest_path_out, UGraph & G )
{
    if (s == d)
    {
        shortest_path_out.push_back(s);
        return;
    }

    std::vector<graph_traits<UGraph>::vertex_descriptor> parent(num_vertices(G));
    typedef graph_traits<UGraph>::vertices_size_type size_type;
    for (size_type p = 0; p < num_vertices(G); ++p)
    {
        parent[p]=p;
    }

    dijkstra_shortest_paths(G,s,predecessor_map(&parent[0]));
    graph_traits<UGraph>::vertex_descriptor next = d;

    do
    {
        shortest_path_out.push_back(next);
        next = parent[next];
    }
    while(next != parent[next]);
    shortest_path_out.push_back(next);
    std::reverse(shortest_path_out.begin(),shortest_path_out.end());
}

bool calculate( UGraph & G, Paths & recalculated_paths )
{
    //all vertices start out as their own parent
    std::vector<graph_traits<UGraph>::vertex_descriptor> parent(num_vertices(G));
    typedef graph_traits<UGraph>::vertices_size_type size_type;
    for (size_type p = 0; p < num_vertices(G); ++p)
    {
        parent[p]=p;
    }

    //traversal of vertex, each as a shortest path source
    typedef graph_traits<UGraph>::vertex_iterator vtx_iter;
    for (std::pair<vtx_iter,vtx_iter> i = vertices(G); i.first != i.second; ++i.first)
    {
        graph_traits<UGraph>::vertex_descriptor s = *i.first;
        if(s == 0)
        {
            continue;
        }

        //compute shortest path
        dijkstra_shortest_paths(G,s,predecessor_map(&parent[0]));

        property_map<UGraph,vertex_index_t>::type vtx_idx_map = get(vertex_index,G);
        for (std::pair<vtx_iter,vtx_iter> j = vertices(G); j.first != j.second; ++j.first)
        {
            graph_traits<UGraph>::vertex_descriptor d = *j.first;
            if (d == 0 || s == d)
            {
                continue;
            }
            if (s != d && s == 0 && d == parent[d])
            {
                return false;
            }
            if (recalculated_paths.find(std::make_pair(s,d)) != recalculated_paths.end())
            {
                continue;
            }
            MSSPF_LOG << s << " -> " << d << " : ";

            std::vector<int> path;
            graph_traits<UGraph>::vertex_descriptor next = d;
            do
            {
                path.push_back(next);
                next = parent[next];
            }
            while(next != parent[next]);
            path.push_back(next);
            std::reverse(path.begin(),path.end());
            recalculated_paths[std::make_pair(s,d)] = path;

            for (size_t i = 0; i < path.size(); ++i)
            {
                MSSPF_LOG <<path[i]<< "->";
            }
            MSSPF_LOG << std::endl;
        }
        MSSPF_LOG << std::endl;
    }
    return true;
}

void Recalculate(UGraph & G, Paths & original_paths,
                 std::vector<graph_traits<UGraph>::edge_descriptor> & failed_edges_v,
                 Paths & rec_paths,
                 Paths & rec_effected_paths
                )
{
    rec_paths = original_paths;
    //find a reroute S->D path
    std::map<std::pair<int,int>,std::pair<int,int>> affected_S_D_paths;
    search_affected(G,original_paths,failed_edges_v,affected_S_D_paths);
    for (auto i = affected_S_D_paths.begin(); i != affected_S_D_paths.end(); ++i)
    {
        int src_router = i->first.first;
        int dst_router = i->first.second;
        int brk_edge_src = i->second.first;
        int brk_dege_dst = i->second.second;

        std::vector<int> shortest_path;
        shortest_path_s_d(src_router,dst_router,shortest_path,G);
        rec_paths[std::make_pair(src_router,dst_router)] = shortest_path;
        rec_effected_paths[std::make_pair(src_router,dst_router)] = shortest_path;
        for (std::size_t i = 0; i != rec_paths[std::make_pair(src_router,dst_router)].size(); ++i)
        {
            MSSPF_LOG << rec_paths[std::make_pair(src_router,dst_router)][i] << "->";
        }
        MSSPF_LOG << std::endl;
    }
}

//MSSPF calculate path
void MSSPF(
    UGraph & G, Paths & original_paths,
    std::vector<graph_traits<UGraph>::edge_descriptor> & failed_edges_v,
    Paths & MSSPF_paths,
    Paths & MSSPF_effected_paths)
{
    MSSPF_paths = original_paths;
    //find a reroute S->D path
    std::map<std::pair<int,int>,std::pair<int,int>> affected_S_D_paths;
    search_affected(G,original_paths,failed_edges_v,affected_S_D_paths);

    for (auto i = affected_S_D_paths.begin(); i != affected_S_D_paths.end(); ++i)
    {
        int src_router = i->first.first;
        int dst_router = i->first.second;
        int brk_edge_src = i->second.first;
        int brk_dege_dst = i->second.second;

        std::vector<int> MSSPF_path,src_brk_src,brk_src_dst;
        shortest_path_s_d(src_router,brk_edge_src,src_brk_src,G);
        shortest_path_s_d(brk_edge_src,dst_router,brk_src_dst,G);

        MSSPF_path = src_brk_src;
        MSSPF_path.insert(MSSPF_path.end(),brk_src_dst.begin() + 1,brk_src_dst.end());
        MSSPF_paths[std::make_pair(src_router,dst_router)] = MSSPF_path;
        MSSPF_effected_paths[std::make_pair(src_router,dst_router)] = MSSPF_path;
        for (std::size_t i = 0; i != MSSPF_paths[std::make_pair(src_router,dst_router)].size(); ++i)
        {
            MSSPF_LOG << MSSPF_paths[std::make_pair(src_router,dst_router)][i] << "->";
        }
        MSSPF_LOG << std::endl;
    }
}

void NotVia( UGraph &G, Paths & original_paths,
             std::vector<graph_traits<UGraph>::edge_descriptor> & failed_edges_v,
             Paths & NotVia_paths ,
             Paths & Notvia_effected_paths)
{
    NotVia_paths = original_paths;

    std::map<std::pair<int,int>,std::pair<int,int>> affected_S_D_paths;
    search_affected(G,original_paths,failed_edges_v,affected_S_D_paths);

    for (auto i = affected_S_D_paths.begin(); i != affected_S_D_paths.end(); ++i)
    {
        int src_router = i->first.first;
        int dst_router = i->first.second;
        int brk_edge_src = i->second.first;
        int brk_dege_dst = i->second.second;
        std::vector<int> not_via_path,src_brk_src,brk_src_dst,brk_dst_dst;
        shortest_path_s_d(src_router,brk_edge_src,src_brk_src,G);
        shortest_path_s_d(brk_edge_src,brk_dege_dst,brk_src_dst,G);
        shortest_path_s_d(brk_dege_dst,dst_router,brk_dst_dst,G);
        if (src_brk_src.size() != 0)
        {
            not_via_path.insert(not_via_path.end(),src_brk_src.begin(),src_brk_src.end());
        }
        if (brk_src_dst.size() != 0)
        {
            not_via_path.insert(not_via_path.end(),brk_src_dst.begin() + 1,brk_src_dst.end());
        }

        if (brk_dst_dst.size() != 0)
        {
            not_via_path.insert(not_via_path.end(),brk_dst_dst.begin() + 1,brk_dst_dst.end());
        }


        NotVia_paths[std::make_pair(src_router,dst_router)] = not_via_path;
        Notvia_effected_paths[std::make_pair(src_router,dst_router)] = not_via_path;
        for (std::size_t i = 0; i != NotVia_paths[std::make_pair(src_router,dst_router)].size(); ++i)
        {
            MSSPF_LOG << NotVia_paths[std::make_pair(src_router,dst_router)][i] << "->";
        }
        MSSPF_LOG << std::endl;
    }
}


void cclt_avg_len(std::string name, Paths & paths )
{
    double length_sum = 0;
    for (auto each_path = paths.begin(); each_path != paths.end(); ++ each_path)
    {
        length_sum += each_path->second.size() - 1;
    }
    std::cout << name << " : " << length_sum/paths.size();
}

void cclt_avg_len(std::string name, std::vector<Paths> & N_paths)
{
    double lth_sum = 0; //total affected path length sum
    double wst_lth_sum = 0;//total path length sum
    double path_ctr = 0; //affected path number
    //traversal each exp in a fixed L size .
    for (auto exp_i = N_paths.begin(); exp_i != N_paths.end(); ++exp_i)
    {
        size_t wst_length = 0;
        path_ctr += exp_i->size();
        for(auto path_i = exp_i->begin(); path_i != exp_i->end(); ++ path_i)
        {
            //traversal each affected path in a exp.
            lth_sum += path_i->second.size() - 1;
            if (path_i->second.size() - 1> wst_length)
            {
                wst_length = path_i->second.size() - 1;
            }
        }
        wst_lth_sum += wst_length;
    }

    std::cout <<"#"<< name <<"\t"<<"avg-l\twrst-l\t" <<std::endl
              <<"\t\t"<<std::setiosflags(std::ios::fixed)<<std::setprecision(3)<<lth_sum/path_ctr<<"\t"<<wst_lth_sum/N_paths.size() << std::endl<<std::endl;
}


void cclt_dstb_path( std::string name, std::vector<Paths> &  N_paths )
{
    std::map<int,int> dstb_map;
    int path_counter = 0;
    for (auto exp_i = N_paths.begin(); exp_i != N_paths.end(); ++exp_i)
    {
        path_counter += exp_i->size();
        for(auto path_i = exp_i->begin(); path_i != exp_i->end(); ++ path_i)
        {
            ++dstb_map[path_i->second.size() - 1];
        }
    }

    std::cout <<"#"<<name <<"-distribution" << std::endl;
    for(auto i = dstb_map.begin(); i != dstb_map.end(); ++i)
    {
        std::cout << i->first << "\t" << i->second <<"\t"<< (double)i->second/path_counter<<std::endl;
    }
    std::cout << std::endl;
}

void compare_length(
    std::vector<Paths> & msspf_N_effected,std::string name_1,
    std::vector<Paths> & notvia_N_effected, std::string name_2 )
{
    std::map<int,int> cmp_map;
    int path_counter = 0;
    for (auto exp_i = 0; exp_i != msspf_N_effected.size(); ++exp_i)
    {
        path_counter += msspf_N_effected.at(exp_i).size();
        for (auto path_i = msspf_N_effected.at(exp_i).begin(); path_i != msspf_N_effected.at(exp_i).end(); ++path_i)
        {
            std::pair<int,int> s_d_pair = path_i->first;
            int path_1_len = path_i->second.size();
            int path_2_len = notvia_N_effected.at(exp_i).at(s_d_pair).size();
            ++cmp_map[path_1_len - path_2_len];
        }
    }

    std::cout <<"#"<< name_1 <<"\t" << name_2 <<std::endl;
    for (auto i = cmp_map.begin(); i != cmp_map.end(); ++i)
    {
        std::cout << i->first << "\t" << i->second <<"\t"<<(double)i->second/path_counter<<std::endl;
    }
    std::cout <<std::endl;
}

typedef struct edge_with_load
{
    Z::Peer_Pair P;
    int load;
    bool operator < (const edge_with_load & other) const
    {
        return load < other.load;
    }
} EDL;

void add_load( Paths & org_paths, std::vector<EDL> & load_ordered )
{
    std::map<Z::peer_pair,int> load_map;
    for (auto i = org_paths.begin(); i != org_paths.end(); ++i)
    {
        for (auto j = i->second.begin(); j != i->second.end() - 1; ++j)
        {
            Z::peer_pair P(*j,*(j+1));
            ++load_map[P];
        }
    }

    for (auto i = load_map.begin(); i != load_map.end(); ++i)
    {
        EDL edl;
        edl.load = i->second;
        edl.P = i->first;
        load_ordered.push_back(edl);
    }
    sort(load_ordered.begin(),load_ordered.end());
}

void add_load(std::vector<Paths> & N_paths, std::map<Z::Peer_Pair,int> & load )
{
    for (auto i = N_paths.begin(); i != N_paths.end(); ++i)
    {
        for (auto j = i->begin(); j != i->end(); ++j)
        {
            for (auto k = j->second.begin(); k != j->second.end() - 1; ++k)
            {
                Z::Peer_Pair P(*k,*(k+1));
                ++load[P];
            }
        }
    }
}

void network_load( Paths & org_paths,
                   std::vector<Paths> & rec_N_paths,
                   std::vector<Paths> & msspf_N_paths,
                   std::vector<Paths> & notvia_N_paths )
{
    std::vector<EDL> org_load_ordered;
    std::map<Z::Peer_Pair,int> rec_load,msspf_load,notvia_load;
    add_load(org_paths,org_load_ordered);
    add_load(rec_N_paths,rec_load);
    add_load(msspf_N_paths,msspf_load);
    add_load(notvia_N_paths,notvia_load);

    std::cout << "#Index	Orgnl	Rclc	MSSPF	Notvia" <<std::endl;
    for (size_t i = 0; i < org_load_ordered.size(); ++i)
    {
        Z::Peer_Pair P = org_load_ordered[i].P;
        std::cout << i <<"\t"
                  << org_load_ordered[i].load << "\t"
                  << (double)rec_load.at(P)/rec_N_paths.size()<< "\t"
                  << (double)msspf_load.at(P)/msspf_N_paths.size()<< "\t"
                  << (double)notvia_load.at(P)/notvia_load.size() << "\t" <<std::endl;

    }
}

int main(int arg,char * argv[])
{
    MSSPF_LOG.open("log.txt");
    UGraph COST239;
    if (import_network_list(argv[1],COST239) != 0)
    {
        return -1;
    }

    MSSPF_LOG << "build graph complete." <<std::endl;

    //calculate original paths
    MSSPF_LOG << "original paths : " <<std::endl;
    Paths org_paths;
    calculate(COST239,org_paths);

    UGraph G = COST239;
    //L is  SRLG size from 1 to 5.
    int L = 5;
    for(int l = 1; l <= L; l++)
    {
        std::cout << "#==================" << std::endl;
        std::cout << "#L = " << l << std::endl;
        MSSPF_LOG << "#=================="<<std::endl;
        MSSPF_LOG << "#L = " << l <<std::endl;
        std::vector<std::vector<graph_traits<UGraph>::edge_descriptor>> failed_edges_v;
        if (l == 1)
        {
            for (auto i = edges(G).first; i != edges(G).second; ++i)
            {
                std::vector<graph_traits<UGraph>::edge_descriptor> v;
                v.push_back(*i);
                if (connected(v,G))
                {
                    failed_edges_v.push_back(v);
                }
            }

        }
        else
        {
            select_failed_paths(l,25,G,failed_edges_v);
        }
        double re_sum = 0, msspf_sum = 0, notvia_sum = 0;
        int re_counter = 0,msspf_counter = 0, notvia_counter = 0;
        std::map<int,int> msspf_VS_notvia;
        std::vector<Paths> rec_N_paths,msspf_N_paths,notvia_N_paths;
        std::vector<Paths> rec_N_effected,msspf_N_effected,notvia_N_effected;

        for (auto i = failed_edges_v.begin(); i != failed_edges_v.end(); ++i)
        {

            //remove the failed paths from original graph: Graph(in),failed edges(in)
            for (auto j = i->begin(); j != i->end(); ++j)
            {
                MSSPF_LOG << "remove edges : " << *j <<std::endl;
            }
            remove_failed_paths(G,*i);

            //recalculate
            MSSPF_LOG << "recalculate : " <<std::endl;
            Paths rec_paths,rec_effected;
            Recalculate(G,org_paths,*i,rec_paths,rec_effected);
            rec_N_paths.push_back(rec_paths);
            rec_N_effected.push_back(rec_effected);


            //MSSPF
            MSSPF_LOG << "#MSSPF" <<std::endl;
            Paths msspf_paths,msspf_effected;
            MSSPF(G,org_paths,*i,msspf_paths,msspf_effected);
            msspf_N_paths.push_back(msspf_paths);
            msspf_N_effected.push_back(msspf_effected);

            //NotVia
            MSSPF_LOG << "#Not Via" <<std::endl;
            Paths notvia_paths,notvia_effected;
            NotVia(G,org_paths,*i,notvia_paths,notvia_effected);
            notvia_N_paths.push_back(notvia_paths);
            notvia_N_effected.push_back(notvia_effected);

            //Recover Graph for next calculation
            G = COST239;
        }

        //statistic

        //average paths length and worst average path length
        std::vector<Paths> org_N_paths;
        org_N_paths.push_back(org_paths);
        cclt_avg_len("Orgnl" , org_N_paths);
        cclt_avg_len("Rcclt" , rec_N_effected);
        cclt_avg_len("MSSPF", msspf_N_effected);
        cclt_avg_len("Notvia", notvia_N_effected);

        //path length distribution
        cclt_dstb_path("Orgnl", org_N_paths);
        cclt_dstb_path("Rcclt", rec_N_paths);
        cclt_dstb_path("MSSPF", msspf_N_paths);
        cclt_dstb_path("Notvia", notvia_N_paths);

        //compare path length
        compare_length(msspf_N_effected,"MSSPF",notvia_N_effected,"Notvia");
        //network load
        network_load(org_paths,rec_N_paths,msspf_N_paths,notvia_N_paths);
        std::cout << std::endl<<std::endl;

    }
    return 0;
}



