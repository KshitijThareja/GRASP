
void graphUtils::build_global_index(int k, int w)
{
    std::cerr << "[M::index] Building Global Index (k=" << k << ", w=" << w << ")..." << std::endl;
    
    // Iterate over all segments
    for (uint32_t i = 0; i < g->n_seg; ++i) {
        gfa_seg_t *s = &g->seg[i];
        if (s->del || s->len < k) continue;
        
        int partition_id = s->partition_id;
        if (partition_id == 0) continue; 
        
        mg128_v mini = {0,0,0};
        // Use mg_sketch for robust minimizer generation
        mg_sketch(0, s->seq, s->len, w, k, i, &mini);
        
        for (size_t j = 0; j < mini.n; ++j) {
            // mini.a[j].x contains the hash in the upper bits? 
            // mg_sketch: p->a[i].x = kMer<<8 | kmerSpan
            // We use the kMer part as the key.
            uint64_t kmer_hash = mini.a[j].x >> 8;
            
            if (global_idx.find(kmer_hash) == global_idx.end()) {
                global_idx[kmer_hash] = partition_id;
            }
        }
        free(mini.a);
    }
    std::cerr << "[M::index] Global Index built. Size: " << global_idx.size() << " minimizers." << std::endl;
}

// Helper to open/close file handles for partitions
std::vector<FILE*> partition_fps;
void graphUtils::open_partition_read_files(const char *prefix)
{
    int max_pid = 0;
    for (auto const& [kmer, pid] : global_idx) {
        if (pid > max_pid) max_pid = pid;
    }
    // Also check segments just in case index is empty?
    // Assuming global_idx covers all partitions.
    
    partition_fps.resize(max_pid + 1, nullptr);
    for (int i = 1; i <= max_pid; ++i) { // Partition 0 is usually "unassigned" or "small", we can write to it too.
        char filename[256];
        sprintf(filename, "%s_%d.fq", prefix, i);
        partition_fps[i] = fopen(filename, "w");
        if (!partition_fps[i]) {
             std::cerr << "[E::dispatch] Failed to open " << filename << std::endl;
        }
    }
    // Handle partition 0
    char filename[256];
    sprintf(filename, "%s_0.fq", prefix);
    partition_fps[0] = fopen(filename, "w");
}

void graphUtils::close_partition_read_files()
{
    for (auto fp : partition_fps) {
        if (fp) fclose(fp);
    }
    partition_fps.clear();
}

int graphUtils::dispatch_read(const char *seq, int len, int k, int w)
{
    if (len < k) return 0;
    
    std::map<int, int> votes;
    
    mg128_v mini = {0,0,0};
    mg_sketch(0, seq, len, w, k, 0, &mini);
    
    for (size_t j = 0; j < mini.n; ++j) {
        uint64_t kmer_hash = mini.a[j].x >> 8;
        if (global_idx.find(kmer_hash) != global_idx.end()) {
            votes[global_idx[kmer_hash]]++;
        }
    }
    free(mini.a);
    
    int best_pid = 0;
    int max_votes = -1;
    
    for (auto const& [pid, count] : votes) {
        if (count > max_votes) {
            max_votes = count;
            best_pid = pid;
        }
    }
    
    return best_pid;
}

void graphUtils::write_read_to_partition(int pid, const char *name, const char *seq, const char *qual)
{
    if (partition_fps.empty()) {
        // Lazy init or error? Let's assume open_partition_read_files was called.
        // If not, just return (or print error once).
        return; 
    }
    if (pid >= 0 && pid < partition_fps.size() && partition_fps[pid]) {
        fprintf(partition_fps[pid], "@%s\n%s\n+\n%s\n", name, seq, qual ? qual : "*"); // Simple FASTQ write
        // If qual is NULL or empty, print *? Or construct fake qual? 
        // mg_bseq1_t has qual.
    }
}
