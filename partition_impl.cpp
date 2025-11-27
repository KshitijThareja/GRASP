
void graphUtils::partition_components(int threshold)
{
    std::cerr << "[M::partition] Starting graph partitioning with threshold " << threshold << std::endl;
    int partitioned_count = 0;

    for (size_t cid = 0; cid < num_cid; cid++)
    {
        int comp_size = conn_comp[cid].size();
        if (comp_size <= threshold) {
            // Assign default partition 0 to all nodes in small components
            for (auto global_id : conn_comp[cid]) {
                g->seg[global_id].partition_id = 0;
            }
            continue;
        }

        partitioned_count++;
        std::cerr << "[M::partition] Partitioning component " << cid << " with " << comp_size << " nodes." << std::endl;

        // Prepare CSR format for METIS
        idx_t nvtxs = comp_size;
        idx_t ncon = 1;
        std::vector<idx_t> xadj(nvtxs + 1);
        std::vector<idx_t> adjncy;
        
        // We need to map local component indices to CSR
        // adj_cc[cid] contains adjacency list with local indices
        
        xadj[0] = 0;
        for (int i = 0; i < nvtxs; ++i) {
            for (auto neighbor : adj_cc[cid][i]) {
                adjncy.push_back(neighbor);
            }
            xadj[i+1] = adjncy.size();
        }

        idx_t nparts = 2; // Start with 2 partitions, or calculate based on size
        if (comp_size > 10000) nparts = 4; // Simple heuristic
        if (comp_size > 100000) nparts = 8;

        std::vector<idx_t> part(nvtxs);
        idx_t objval;
        
        // METIS options
        idx_t options[METIS_NOPTIONS];
        METIS_SetDefaultOptions(options);
        options[METIS_OPTION_SEED] = 42;

        int ret = METIS_PartGraphKway(&nvtxs, &ncon, xadj.data(), adjncy.data(), 
                                      NULL, NULL, NULL, &nparts, NULL, NULL, options, &objval, part.data());

        if (ret != METIS_OK) {
            std::cerr << "[E::partition] METIS failed for component " << cid << std::endl;
            // Fallback: assign 0
            for (auto global_id : conn_comp[cid]) {
                g->seg[global_id].partition_id = 0;
            }
        } else {
            // Assign partition IDs back to graph segments
            for (int i = 0; i < nvtxs; ++i) {
                int global_id = component_idx[cid][i];
                g->seg[global_id].partition_id = part[i];
            }
            std::cerr << "[M::partition] Component " << cid << " partitioned into " << nparts << " parts. Cut: " << objval << std::endl;
        }
    }
    std::cerr << "[M::partition] Partitioning complete. Partitioned " << partitioned_count << " components." << std::endl;
}
