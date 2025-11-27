
void graphUtils::identify_halo_nodes()
{
    std::cerr << "[M::halo] Identifying halo nodes..." << std::endl;
    int halo_count = 0;
    
    // Iterate over all segments
    for (uint32_t v = 0; v < n_vtx; ++v) {
        int seg_id = v >> 1;
        int my_part = g->seg[seg_id].partition_id;
        
        int n_edges = gfa_arc_n(g, v);
        gfa_arc_t *av = gfa_arc_a(g, v);
        
        for (int i = 0; i < n_edges; ++i) {
            uint32_t w = av[i].w;
            int neighbor_seg_id = w >> 1;
            int neighbor_part = g->seg[neighbor_seg_id].partition_id;
            
            if (my_part != neighbor_part) {
                // Boundary detected!
                if (!g->seg[seg_id].is_halo) {
                    g->seg[seg_id].is_halo = 1;
                    halo_count++;
                }
                if (!g->seg[neighbor_seg_id].is_halo) {
                    g->seg[neighbor_seg_id].is_halo = 1;
                    halo_count++;
                }
            }
        }
    }
    std::cerr << "[M::halo] Identified " << halo_count << " halo nodes." << std::endl;
}

void graphUtils::generate_partition_files(const char *prefix)
{
    // Find max partition id
    int max_pid = 0;
    for (uint32_t i = 0; i < g->n_seg; ++i) {
        if (g->seg[i].partition_id > max_pid) {
            max_pid = g->seg[i].partition_id;
        }
    }
    
    std::cerr << "[M::output] Generating " << max_pid + 1 << " partition files..." << std::endl;
    
    for (int pid = 0; pid <= max_pid; ++pid) {
        char filename[256];
        sprintf(filename, "%s_%d.gfa", prefix, pid);
        FILE *fp = fopen(filename, "w");
        if (!fp) {
            std::cerr << "[E::output] Cannot open file " << filename << " for writing." << std::endl;
            continue;
        }
        
        gfa_print_partition(g, fp, pid);
        
        fclose(fp);
    }
    std::cerr << "[M::output] Partition files generated." << std::endl;
}
