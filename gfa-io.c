#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include "kstring.h"
#include "gfa-priv.h"

#include "kseq.h"
KSTREAM_INIT(gzFile, gzread, 65536)

/***********
 * Tag I/O *
 ***********/

int gfa_aux_parse(char *s, uint8_t **data, int *max)
{
	char *q, *p;
	kstring_t str;
	if (s == 0) return 0;
	str.l = 0, str.m = *max, str.s = (char*)*data;
	if (*s == '\t') ++s;
	for (p = q = s;; ++p) {
		if (*p == 0 || *p == '\t') {
			int c = *p;
			*p = 0;
			if (p - q >= 5 && q[2] == ':' && q[4] == ':' && (q[3] == 'A' || q[3] == 'i' || q[3] == 'f' || q[3] == 'Z' || q[3] == 'B')) {
				int type = q[3];
				kputsn_(q, 2, &str);
				q += 5;
				if (type == 'A') {
					kputc_('A', &str);
					kputc_(*q, &str);
				} else if (type == 'i') {
					int32_t x;
					x = strtol(q, &q, 10);
					kputc_(type, &str); kputsn_((char*)&x, 4, &str);
				} else if (type == 'f') {
					float x;
					x = strtod(q, &q);
					kputc_('f', &str); kputsn_(&x, 4, &str);
				} else if (type == 'Z') {
					kputc_('Z', &str); kputsn_(q, p - q + 1, &str); // note that this include the trailing NULL
				} else if (type == 'B') {
					type = *q++; // q points to the first ',' following the typing byte
					if (p - q >= 2 && (type == 'c' || type == 'C' || type == 's' || type == 'S' || type == 'i' || type == 'I' || type != 'f')) {
						int32_t n;
						char *r;
						for (r = q, n = 0; *r; ++r)
							if (*r == ',') ++n;
						kputc_('B', &str); kputc_(type, &str); kputsn_(&n, 4, &str);
						// TODO: to evaluate which is faster: a) aligned array and then memmove(); b) unaligned array; c) kputsn_()
						if (type == 'c')      while (q + 1 < p) { int8_t   x = strtol(q + 1, &q, 0); kputc_(x, &str); }
						else if (type == 'C') while (q + 1 < p) { uint8_t  x = strtol(q + 1, &q, 0); kputc_(x, &str); }
						else if (type == 's') while (q + 1 < p) { int16_t  x = strtol(q + 1, &q, 0); kputsn_(&x, 2, &str); }
						else if (type == 'S') while (q + 1 < p) { uint16_t x = strtol(q + 1, &q, 0); kputsn_(&x, 2, &str); }
						else if (type == 'i') while (q + 1 < p) { int32_t  x = strtol(q + 1, &q, 0); kputsn_(&x, 4, &str); }
						else if (type == 'I') while (q + 1 < p) { uint32_t x = strtol(q + 1, &q, 0); kputsn_(&x, 4, &str); }
						else if (type == 'f') while (q + 1 < p) { float    x = strtod(q + 1, &q);    kputsn_(&x, 4, &str); }
					}
				} // should not be here, as we have tested all types
			}
			q = p + 1;
			if (c == 0) break;
		}
	}
	if (str.l > 0 && str.l == str.m) ks_resize(&str, str.l + 1);
	if (str.s) str.s[str.l] = 0;
	*max = str.m, *data = (uint8_t*)str.s;
	return str.l;
}

int gfa_aux_format(int l_aux, const uint8_t *aux, char **t, int *max)
{
	kstring_t str;
	const uint8_t *s = aux;
	str.l = 0, str.s = *t, str.m = *max;
	while (s < aux + l_aux) {
		uint8_t type, key[2];
		key[0] = s[0]; key[1] = s[1];
		s += 2; type = *s++;
		kputc('\t', &str); kputsn((char*)key, 2, &str); kputc(':', &str);
		if (type == 'A') { kputsn("A:", 2, &str); kputc(*s, &str); ++s; }
		else if (type == 'i') { kputsn("i:", 2, &str); kputw(*(int32_t*)s, &str); s += 4; }
		else if (type == 'f') { ksprintf(&str, "f:%g", *(float*)s); s += 4; }
		else if (type == 'Z') { kputc(type, &str); kputc(':', &str); while (*s) kputc(*s++, &str); ++s; }
		else if (type == 'B') {
			uint8_t sub_type = *(s++);
			int32_t i, n;
			memcpy(&n, s, 4);
			s += 4; // no point to the start of the array
			kputsn("B:", 2, &str); kputc(sub_type, &str); // write the typing
			for (i = 0; i < n; ++i) { // FIXME: for better performance, put the loop after "if"
				kputc(',', &str);
				if ('c' == sub_type)      { kputw(*(int8_t*)s, &str); ++s; }
				else if ('C' == sub_type) { kputw(*(uint8_t*)s, &str); ++s; }
				else if ('s' == sub_type) { kputw(*(int16_t*)s, &str); s += 2; }
				else if ('S' == sub_type) { kputw(*(uint16_t*)s, &str); s += 2; }
				else if ('i' == sub_type) { kputw(*(int32_t*)s, &str); s += 4; }
				else if ('I' == sub_type) { kputuw(*(uint32_t*)s, &str); s += 4; }
				else if ('f' == sub_type) { ksprintf(&str, "%g", *(float*)s); s += 4; }
			}
		}
	}
	*t = str.s, *max = str.m;
	return str.l;
}

/****************
 * Line parsers *
 ****************/

int gfa_parse_S(gfa_t *g, char *s)
{
	int i, is_ok = 0;
	char *p, *q, *seg = 0, *seq = 0, *rest = 0;
	uint32_t sid, len = 0;
	for (i = 0, p = q = s + 2;; ++p) {
		if (*p == 0 || *p == '\t') {
			int c = *p;
			*p = 0;
			if (i == 0) seg = q;
			else if (i == 1) {
				seq = q[0] == '*'? 0 : gfa_strdup(q);
				is_ok = 1, rest = c? p + 1 : 0;
				break;
			}
			++i, q = p + 1;
			if (c == 0) break;
		}
	}
	if (is_ok) { // all mandatory fields read
		int l_aux, m_aux = 0, LN = -1;
		uint8_t *aux = 0, *s_LN = 0;
		gfa_seg_t *s;
		l_aux = gfa_aux_parse(rest, &aux, &m_aux); // parse optional tags
		s_LN = l_aux? gfa_aux_get(l_aux, aux, "LN") : 0;
		if (s_LN && s_LN[0] == 'i') {
			LN = *(int32_t*)(s_LN + 1);
			l_aux = gfa_aux_del(l_aux, aux, s_LN);
		}
		if (seq == 0) {
			if (LN >= 0) len = LN;
		} else len = strlen(seq);
		if (LN >= 0 && len != LN && gfa_verbose >= 2)
			fprintf(stderr, "[W] for segment '%s', LN:i:%d tag is different from sequence length %d\n", seg, LN, len);
		sid = gfa_add_seg(g, seg);
		s = &g->seg[sid];
		s->len = len, s->seq = seq;
		if (l_aux > 0) {
			uint8_t *s_SN = 0, *s_SO = 0, *s_SR = 0;
			s_SN = gfa_aux_get(l_aux, aux, "SN");
			if (s_SN && *s_SN == 'Z') { // then parse stable tags
				s->snid = gfa_sseq_add(g, (char*)(s_SN + 1)), s->soff = 0;
				l_aux = gfa_aux_del(l_aux, aux, s_SN);
				s_SO = gfa_aux_get(l_aux, aux, "SO");
				if (s_SO && *s_SO == 'i') {
					s->soff = *(int32_t*)(s_SO + 1);
					l_aux = gfa_aux_del(l_aux, aux, s_SO);
				}
			}
			s_SR = gfa_aux_get(l_aux, aux, "SR");
			if (s_SR && *s_SR == 'i') {
				s->rank = *(int32_t*)(s_SR + 1);
				if (s->rank > g->max_rank) g->max_rank = s->rank;
				l_aux = gfa_aux_del(l_aux, aux, s_SR);
			}
			gfa_sseq_update(g, s);
		}
		if (l_aux > 0)
			s->aux.m_aux = m_aux, s->aux.l_aux = l_aux, s->aux.aux = aux;
		else if (aux)
			free(aux);
	} else return -1;
	return 0;
}

int gfa_parse_L(gfa_t *g, char *s)
{
	int i, oriv = -1, oriw = -1, is_ok = 0;
	char *p, *q, *segv = 0, *segw = 0, *rest = 0;
	int32_t ov = INT32_MAX, ow = INT32_MAX;
	for (i = 0, p = q = s + 2;; ++p) {
		if (*p == 0 || *p == '\t') {
			int c = *p;
			*p = 0;
			if (i == 0) {
				segv = q;
			} else if (i == 1) {
				if (*q != '+' && *q != '-') return -2;
				oriv = (*q != '+');
			} else if (i == 2) {
				segw = q;
			} else if (i == 3) {
				if (*q != '+' && *q != '-') return -2;
				oriw = (*q != '+');
			} else if (i == 4) {
				if (*q == '*') {
					ov = ow = 0;
				} else if (*q == ':') {
					ov = INT32_MAX;
					ow = isdigit(*(q+1))? strtol(q+1, &q, 10) : INT32_MAX;
				} else if (isdigit(*q)) {
					char *r;
					ov = strtol(q, &r, 10);
					if (isupper(*r)) { // CIGAR
						ov = ow = 0;
						do {
							long l;
							l = strtol(q, &q, 10);
							if (*q == 'M' || *q == 'D' || *q == 'N') ov += l;
							if (*q == 'M' || *q == 'I' || *q == 'S') ow += l;
							++q;
						} while (isdigit(*q));
					} else if (*r == ':') { // overlap lengths
						ow = isdigit(*(r+1))? strtol(r+1, &r, 10) : INT32_MAX;
					} else break;
				} else break;
				is_ok = 1, rest = c? p + 1 : 0;
				break;
			}
			++i, q = p + 1;
			if (c == 0) break;
		}
	}
	if (i == 4 && is_ok == 0) ov = ow = 0, is_ok = 1; // no overlap field
	if (is_ok) {
		uint32_t v, w;
		int l_aux, m_aux = 0;
		uint8_t *aux = 0;
		gfa_arc_t *arc;
		v = gfa_add_seg(g, segv) << 1 | oriv;
		w = gfa_add_seg(g, segw) << 1 | oriw;
		arc = gfa_add_arc1(g, v, w, ov, ow, -1, 0);
		l_aux = gfa_aux_parse(rest, &aux, &m_aux); // parse optional tags
		if (l_aux) {
			gfa_aux_t *a = &g->link_aux[arc->link_id];
			uint8_t *s_L1, *s_L2, *s_SR;
			a->l_aux = l_aux, a->m_aux = m_aux, a->aux = aux;
			s_SR = gfa_aux_get(a->l_aux, a->aux, "SR");
			if (s_SR && s_SR[0] == 'i') {
				arc->rank = *(int32_t*)(s_SR+1);
				a->l_aux = gfa_aux_del(a->l_aux, a->aux, s_SR);
			}
			s_L1 = gfa_aux_get(a->l_aux, a->aux, "L1");
			if (s_L1) {
				if (ov != INT32_MAX && s_L1[0] == 'i')
					g->seg[v>>1].len = g->seg[v>>1].len > ov + *(int32_t*)(s_L1+1)? g->seg[v>>1].len : ov + *(int32_t*)(s_L1+1);
				a->l_aux = gfa_aux_del(a->l_aux, a->aux, s_L1);
			}
			s_L2 = gfa_aux_get(a->l_aux, a->aux, "L2");
			if (s_L2) {
				if (ow != INT32_MAX && s_L2[0] == 'i')
					g->seg[w>>1].len = g->seg[w>>1].len > ow + *(int32_t*)(s_L2+1)? g->seg[w>>1].len : ow + *(int32_t*)(s_L2+1);
				a->l_aux = gfa_aux_del(a->l_aux, a->aux, s_L2);
			}
			if (a->l_aux == 0) {
				free(a->aux);
				a->aux = 0, a->m_aux = 0;
			}
		}
	} else return -1;
	return 0;
}

static gfa_seg_t *gfa_parse_fa_hdr(gfa_t *g, char *s)
{
	int32_t i;
	char buf[16];
	gfa_seg_t *seg;
	for (i = 0; s[i]; ++i)
		if (isspace(s[i])) break;
	s[i] = 0;
	sprintf(buf, "s%d", g->n_seg + 1);
	i = gfa_add_seg(g, buf);
	seg = &g->seg[i];
	seg->snid = gfa_sseq_add(g, s + 1);
	seg->soff = seg->rank = 0;
	return seg;
}

static void gfa_update_fa_seq(gfa_t *g, gfa_seg_t *seg, int32_t l_seq, const char *seq)
{
	if (seg == 0) return;
	seg->seq = gfa_strdup(seq);
	seg->len = l_seq;
	gfa_sseq_update(g, seg);
}

/****************
 * User-end I/O *
 ****************/

gfa_t *gfa_read(const char *fn)
{
	gzFile fp;
	gfa_t *g;
	kstring_t s = {0,0,0}, fa_seq = {0,0,0};
	kstream_t *ks;
	int dret, is_fa = 0;
	gfa_seg_t *fa_seg = 0;
	uint64_t lineno = 0;

	fp = fn && strcmp(fn, "-")? gzopen(fn, "r") : gzdopen(0, "r");
	if (fp == 0) return 0;
	ks = ks_init(fp);
	g = gfa_init();
	while (ks_getuntil(ks, KS_SEP_LINE, &s, &dret) >= 0) {
		int ret = 0;
		++lineno;
		if (s.l > 0 && s.s[0] == '>') { // FASTA header
			is_fa = 1;
			if (fa_seg) gfa_update_fa_seq(g, fa_seg, fa_seq.l, fa_seq.s);
			fa_seg = gfa_parse_fa_hdr(g, s.s);
			fa_seq.l = 0;
		} else if (is_fa) { // FASTA mode
			if (s.l >= 3 && s.s[1] == '\t') { // likely a GFA line
				gfa_update_fa_seq(g, fa_seg, fa_seq.l, fa_seq.s); // finalize fa_seg
				fa_seg = 0;
				is_fa = 0;
			} else kputsn(s.s, s.l, &fa_seq); // likely a FASTA sequence line
		}
		if (is_fa) continue;
		if (s.l < 3 || s.s[1] != '\t') continue; // empty line
		if (s.s[0] == 'S') ret = gfa_parse_S(g, s.s);
		else if (s.s[0] == 'L') ret = gfa_parse_L(g, s.s);
		if (ret < 0 && gfa_verbose >= 1)
			fprintf(stderr, "[E] invalid %c-line at line %ld (error code %d)\n", s.s[0], (long)lineno, ret);
	}
	if (is_fa && fa_seg) gfa_update_fa_seq(g, fa_seg, fa_seq.l, fa_seq.s);
	free(fa_seq.s);
	free(s.s);
	gfa_finalize(g);
	ks_destroy(ks);
	gzclose(fp);
	return g;
}

void gfa_print(const gfa_t *g, FILE *fp, int flag)
{
	uint32_t i;
	uint64_t k;
	for (i = 0; i < g->n_seg; ++i) {
		const gfa_seg_t *s = &g->seg[i];
		if (s->del) continue;
		fprintf(fp, "S\t%s\t", s->name);
		if (s->seq && !(flag & GFA_O_NO_SEQ)) fputs(s->seq, fp);
		else fputc('*', fp);
		fprintf(fp, "\tLN:i:%d", s->len);
		if (s->snid >= 0 && s->soff >= 0)
			fprintf(fp, "\tSN:Z:%s\tSO:i:%d", g->sseq[s->snid].name, s->soff);
		if (s->rank >= 0)
			fprintf(fp, "\tSR:i:%d", s->rank);
		if (s->utg && s->utg->n) fprintf(fp, "\tRC:i:%d\tlc:i:%d", s->utg->n, s->utg->len_comp);
		if (s->aux.l_aux > 0) {
			char *t = 0;
			int max = 0;
			gfa_aux_format(s->aux.l_aux, s->aux.aux, &t, &max);
			fputs(t, fp);
			free(t);
		}
		fputc('\n', fp);
		if (s->utg && s->utg->n) {
			uint32_t j, l;
			for (j = l = 0; j < s->utg->n; ++j) {
				const gfa_utg_t *u = s->utg;
				fprintf(fp, "A\t%s\t%d\t%c\t%s\t%d\t%d\n", s->name, l, "+-"[u->a[j]>>32&1], u->name[j], (int32_t)(u->r[j]>>32), (int32_t)u->r[j]);
				l += (uint32_t)u->a[j];
			}
		}
	}
	for (k = 0; k < g->n_arc; ++k) {
		const gfa_arc_t *a = &g->arc[k];
		const gfa_aux_t *aux = a->link_id < g->n_arc? &g->link_aux[a->link_id] : 0;
		if (a->del || a->comp) continue;
		fprintf(fp, "L\t%s\t%c\t%s\t%c", g->seg[a->v_lv>>33].name, "+-"[a->v_lv>>32&1], g->seg[a->w>>1].name, "+-"[a->w&1]);
		if (!(flag & GFA_O_OV_EXT)) {
			fprintf(fp, "\t%dM", a->ov < a->ow? a->ov : a->ow);
		} else {
			if (a->ov == a->ow) fprintf(fp, "\t%dM", a->ov);
			else fprintf(fp, "\t%d:%d", a->ov, a->ow);
		}
		if (a->rank >= 0) fprintf(fp, "\tSR:i:%d", a->rank);
		fprintf(fp, "\tL1:i:%d", gfa_arc_len(*a));
		fprintf(fp, "\tL2:i:%d", gfa_arc_lw(g, *a));
		if (aux && aux->l_aux) {
			char *t = 0;
			int max = 0;
			gfa_aux_format(aux->l_aux, aux->aux, &t, &max);
			if (t) fputs(t, fp);
			free(t);
		}
		fputc('\n', fp);
	}
}

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
