
void gfa_print_partition(const gfa_t *g, FILE *fp, int partition_id)
{
    uint32_t i;
    for (i = 0; i < g->n_seg; ++i) {
        const gfa_seg_t *s = &g->seg[i];
        if (s->del) continue;
        
        // Check if segment belongs to this partition (Core) or is a Halo node connected to it
        int include = 0;
        if (s->partition_id == partition_id) {
            include = 1;
        } else if (s->is_halo) {
            // Check neighbors to see if any are in this partition
            uint32_t k;
            for (k = 0; k < 2; ++k) {
                uint32_t v = i<<1 | k;
                uint32_t nv = gfa_arc_n(g, v);
                gfa_arc_t *av = gfa_arc_a(g, v);
                for (uint32_t j = 0; j < nv; ++j) {
                    uint32_t w = av[j].w;
                    if (g->seg[w>>1].partition_id == partition_id) {
                        include = 1;
                        break;
                    }
                }
                if (include) break;
            }
        }
        
        if (!include) continue;

        // Print Segment (S-line)
        fprintf(fp, "S\t%s\t", s->name);
        if (s->seq && !(0 & GFA_O_NO_SEQ)) fputs(s->seq, fp);
        else fputc('*', fp);
        fprintf(fp, "\tLN:i:%d", s->len);
        if (s->snid >= 0 && s->soff >= 0)
            fprintf(fp, "\tSN:Z:%s\tSO:i:%d", g->sseq[s->snid].name, s->soff);
        if (s->rank >= 0)
            fprintf(fp, "\tSR:i:%d", s->rank);
        // Add partition tag
        fprintf(fp, "\tpt:i:%d", s->partition_id);
        if (s->is_halo) fprintf(fp, "\thl:i:1");
        
        if (s->aux.l_aux > 0) {
            char *t = 0;
            int max = 0;
            gfa_aux_format(s->aux.l_aux, s->aux.aux, &t, &max);
            if (t) fputs(t, fp);
            free(t);
        }
        fputc('\n', fp);

        // Print Arcs (L-lines) originating from this segment
        // Only if the target is also included in this partition view
        uint32_t k;
        for (k = 0; k < 2; ++k) {
            uint32_t v = i<<1 | k;
            uint32_t nv = gfa_arc_n(g, v);
            gfa_arc_t *av = gfa_arc_a(g, v);
            for (uint32_t j = 0; j < nv; ++j) {
                gfa_arc_t *a = &av[j];
                if (a->del) continue;
                
                // Check if target w is included
                uint32_t w_seg_id = a->w >> 1;
                int include_w = 0;
                if (g->seg[w_seg_id].partition_id == partition_id) {
                    include_w = 1;
                } else if (g->seg[w_seg_id].is_halo) {
                     // Check if w is connected to partition_id (it is, via v!)
                     // But we only want to print edges that are relevant to this partition.
                     // If v is in partition, and w is halo, we include the edge.
                     // If v is halo, and w is in partition, we include the edge.
                     // If both are halo, we might include it if they are both relevant? 
                     // For simplicity: include edge if BOTH endpoints are "included" in this partition view.
                     
                     // We already know v is included. Now check w.
                     // Re-run the include logic for w? Or simpler:
                     // Just check if w is in partition OR w is halo connected to partition.
                     // Since v is in partition (or halo connected), and v->w exists...
                     // If v is CORE, w is definitely relevant (CORE or HALO).
                     // If v is HALO, w MUST be CORE for this edge to be relevant? 
                     // No, if v is HALO and w is HALO, they might be part of the boundary.
                     // Let's stick to: Include edge if both u and v are in the set of nodes for this partition.
                     
                     // We need to check if w would be selected by the node selection logic above.
                     // This is slightly expensive to re-compute.
                     // Optimization: Pre-compute a bitmask or flag for "in_current_partition"? 
                     // For now, just re-check basic condition:
                     if (g->seg[w_seg_id].partition_id == partition_id) include_w = 1;
                     else {
                         // w is halo. Is it connected to ANY core node of this partition?
                         // Yes, if v is core.
                         if (s->partition_id == partition_id) include_w = 1;
                         else {
                             // v is halo. w is halo. 
                             // Are they connected to the core?
                             // We only include the subgraph induced by (Core U Halo).
                             // So yes, if w is halo, it is in the set.
                             // Wait, "Halo" means "connected to a different partition".
                             // A node is in the set if (PartID == PID) OR (IsHalo AND ConnectedTo(PID)).
                             // We know v is in the set.
                             // Is w in the set?
                             // We need to check if w is connected to PID.
                             // We know w is connected to v.
                             // If v is Core(PID), then w is connected to PID. So w is in set.
                             // If v is Halo, and w is Halo... w might NOT be connected to PID directly?
                             // But we want the induced subgraph. So we must check if w is valid for this partition.
                             
                             // Let's duplicate the check logic for w
                             uint32_t wk;
                             for (wk = 0; wk < 2; ++wk) {
                                uint32_t wv = w_seg_id<<1 | wk;
                                uint32_t wnv = gfa_arc_n(g, wv);
                                gfa_arc_t *wav = gfa_arc_a(g, wv);
                                for (uint32_t wj = 0; wj < wnv; ++wj) {
                                    if (g->seg[wav[wj].w>>1].partition_id == partition_id) {
                                        include_w = 1;
                                        break;
                                    }
                                }
                                if (include_w) break;
                             }
                         }
                     }
                }
                
                if (!include_w) continue;

                // Print Arc
                fprintf(fp, "L\t%s\t%c\t%s\t%c", g->seg[i].name, "+-"[v&1], g->seg[w_seg_id].name, "+-"[a->w&1]);
                if (!(0 & GFA_O_OV_EXT)) {
                    fprintf(fp, "\t%dM", a->ov < a->ow? a->ov : a->ow);
                } else {
                    if (a->ov == a->ow) fprintf(fp, "\t%dM", a->ov);
                    else fprintf(fp, "\t%d:%d", a->ov, a->ow);
                }
                if (a->rank >= 0) fprintf(fp, "\tSR:i:%d", a->rank);
                fprintf(fp, "\tL1:i:%d", gfa_arc_len(*a));
                fprintf(fp, "\tL2:i:%d", gfa_arc_lw(g, *a));
                fputc('\n', fp);
            }
        }
    }
}
