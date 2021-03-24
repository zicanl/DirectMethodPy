import math

EXP = 2.71828182845904523536028747135
MSKSEG = [0x0000000000000000,
          0xFF00000000000000,
          0x00FF000000000000,
          0x0000FF0000000000,
          0x000000FF00000000,
          0x00000000FF000000,
          0x0000000000FF0000,
          0x000000000000FF00,
          0x00000000000000FF]

DELSEG = [0x0000000000000000,
          0x0000FFFFFFFFFFFF,
          0xFF0000FFFFFFFFFF,
          0xFFFF0000FFFFFFFF,
          0xFFFFFF0000FFFFFF,
          0xFFFFFFFF0000FFFF,
          0xFFFFFFFFFF0000FF,
          0xFFFFFFFFFFFF0000]

MSKBIT = [0x0000000000000000,
          0x8080808080808080,
          0x4040404040404040,
          0x2020202020202020,
          0x1010101010101010,
          0x0808080808080808,
          0x0404040404040404,
          0x0202020202020202,
          0x0101010101010101]

DELBIT = [0x0000000000000000,
          0x3F3F3F3F3F3F3F3F,
          0x9F9F9F9F9F9F9F9F,
          0xCFCFCFCFCFCFCFCF,
          0xE7E7E7E7E7E7E7E7,
          0xF3F3F3F3F3F3F3F3,
          0xF9F9F9F9F9F9F9F9,
          0xFCFCFCFCFCFCFCFC]


class Bender:
    def __init__(self, n, r, c, rand):
        if c is not None:
            self.directed = True
            self.n = n
            self.r = r
            self.c = c
            tmp = 0
            self.rand = rand
            self.f = 0
            self.sqsum_r = 0
            self.sqsum_c = 0
            self.logfac_r = [0 for _ in range(n)]
            self.logfac_c = [0 for _ in range(n)]
            self.maxdeg = 0
            for i in range(n):
                tmp += (c[i]) * (r[i])
                self.logfac_r[i] = math.lgamma(r[i] + 1)
                self.logfac_c[i] = math.lgamma(c[i] + 1)
                if r[i] > self.maxdeg:
                    self.maxdeg = r[i]
                if c[i] > self.maxdeg:
                    self.maxdeg = c[i]
                self.f += r[i]
                self.sqsum_r += r[i] * r[i] - r[i]
                self.sqsum_c += c[i] * c[i] - c[i]

            self.logfac_f = [0 for i in range(64)]
            for i in range(64):
                self.logfac_f[i] = math.lgamma(self.f - i + 1)
            self.logfac_cache = [0 for i in range(self.maxdeg + 1)]
            for i in range(self.maxdeg + 1):
                self.logfac_cache[i] = math.lgamma(i + 1)
            self.b = tmp / self.f
        else:
            self.directed = False
            self.n = n
            self.c = None
            self.r = r
            self.b = 0
            self.rand = rand
            self.f = 0
            self.sqsum_r = 0
            self.logfac_r = [0 for _ in range(n)]
            for i in range(n):
                self.logfac_r[i] = math.lgamma(r[i] + 1)
                self.f += r[i]
                self.sqsum_r += r[i] * r[i] - r[i]

    def motif_sample_log(self, motif, k, num_samples):
        test = []
        results = []
        result_buffer = []
        BUF_SIZE = 1000000
        odg = []
        idg = []
        self._noautoms(odg, idg, motif, k, self.directed, results)
        num_graphs = len(odg) // k
        odeg = [0 for _ in range(k)]
        ideg = [0 for _ in range(k)]
        pos = [0 for _ in range(k)]

        candidate_indices = []
        for i in range(self.n):
            if self.directed:
                for j in range(k):
                    if self.r[i] >= odg[j] and self.c[i] >= idg[j]:
                        candidate_indices.append(i)
                        break
            else:
                for j in range(k):
                    if self.r[i] >= odg[j]:
                        candidate_indices.append(i)
                        break

        new_n = len(candidate_indices)
        taken = [False for i in range(new_n)]

        maximum = 0.0
        cand = 0
        if new_n >= k:
            for i in range(num_samples):
                for j in range(k):
                    cand = int(self.rand() * 10 ** 8) % new_n
                    while taken[cand]:
                        cand = int(self.rand() * 10 ** 8) % new_n
                    taken[cand] = True
                    pos[j] = cand

                for j in range(k):
                    taken[pos[j]] = False
                    pos[j] = candidate_indices[pos[j]]

                for j in range(num_graphs):
                    if self.directed:
                        for l in range(k):
                            odeg[l] = odg[j * k + l]
                            ideg[l] = idg[j * k + l]
                        self._calc_motif_dir(pos, odeg, ideg, k, results)
                    else:
                        for l in range(k):
                            odeg[l] = odg[j * k + l]
                        self._calc_motif_undir(pos, odeg, k, results, test)

                    if len(results) > BUF_SIZE:
                        sum_ = 0.0
                        l = 0
                        results.sort()
                        if len(result_buffer) == 0:
                            maximum = results[len(results) - 1]
                        while l != BUF_SIZE:
                            sum_ += math.exp(results[l] - maximum)
                            l += 1
                        results = results[BUF_SIZE + 1:]
                        result_buffer.append(sum_)

        p = 1.0
        for i in range(k):
            p *= (new_n - i) / (self.n - i)

        results.sort()
        ret = 0.0
        if len(results) > 0:
            if len(result_buffer) == 0:
                maximum = results[-1]
            sum_ = 0
            for i in range(len(results)):
                sum_ += math.exp(results[i] - maximum)
            for i in range(len(result_buffer)):
                sum_ += result_buffer[i]
            ret = maximum + math.log(sum_) - math.log(num_samples / p)
        else:
            # ret = math.sqrt(-1.0)
            ret = -float('nan')
        return ret

    def _getcanonical(self, g, k):
        gr = g
        tmp1, tmp2 = 0, 0
        c = [0 for i in range(k + 1)]
        o = [0 for i in range(k + 1)]
        j = k
        s, q = 0, 0
        t1, t2 = 0, 0
        canon = g
        for i in range(k + 1):
            c[i] = 0
            o[i] = 1
        while j != 0:
            q = c[j] + o[j]
            if q < 0:
                o[j] = -o[j]
                j -= 1
            elif q == j:
                s += 1
                o[j] = -o[j]
                j -= 1
            elif q > -1:
                t1 = j - q + s
                t2 = j - c[j] + s
                if t1 > t2:
                    t1 ^= t2
                    t2 ^= t1
                    t1 ^= t2
                tmp1 = gr & MSKSEG[t1]
                tmp2 = gr & MSKSEG[t2]
                gr &= DELSEG[t1]
                gr |= (tmp1 >> 8)
                gr |= (tmp2 << 8)
                tmp1 = gr & MSKBIT[t1]
                tmp2 = gr & MSKBIT[t2]
                gr &= DELBIT[t1]
                gr |= (tmp1 >> 1)
                gr |= (tmp2 << 1)
                if gr > canon:
                    canon = gr
                c[j] = q
                j = k
                s = 0
        return canon

    def _noautoms(self, odg, idg, g, k, directed, result):
        gr = int(0)
        already = set()
        odeg = [0 for i in range(k)]
        ideg = [0 for i in range(k)]
        for i in range(k):
            for j in range(k):
                if (g & (int(1) << i + j * k)) > 0:
                    gr = SET(gr, j, i)
                    ideg[i] += 1
                    odeg[j] += 1
        tmp1, tmp2 = 0, 0
        c = [0 for i in range(k+1)]
        o = [1 for i in range(k+1)]
        j = k
        s = 0
        q, t1, t2 = 0, 0, 0
        already.add(gr)
        for i in range(k):
            idg.append(ideg[i])
            odg.append(odeg[i])

        while j != 0:
            q = c[j] + o[j]
            if q < 0:
                o[j] = -o[j]
                j -= 1
            elif q == j:
                s += 1
                o[j] = -o[j]
                j -= 1
            elif q > -1:
                t1 = (j-q+s)
                t2 = (j-c[j]+s)
                if t1 > t2:
                    t1 ^= t2
                    t2 ^= t1
                    t1 ^= t2
                tmp1 = ideg[t1 - 1]
                ideg[t1 - 1] = ideg[t2 - 1]
                ideg[t2 - 1] = tmp1
                tmp1 = odeg[t1 - 1]
                odeg[t1 - 1] = odeg[t2 - 1]
                odeg[t2 - 1] = tmp1
                tmp1 = gr & MSKSEG[t1]
                tmp2 = gr & MSKSEG[t2]
                gr &= DELSEG[t1]
                gr |= (tmp1 >> 8)
                gr |= (tmp2 << 8)
                tmp1 = gr & MSKBIT[t1]
                tmp2 = gr & MSKBIT[t2]
                gr &= DELBIT[t1]
                gr |= (tmp1 >> 1)
                gr |= (tmp2 << 1)
                if gr not in already:
                    already.add(gr)
                    for i in range(k):
                        idg.append(ideg[i])
                        odg.append(odeg[i])
                c[j] = q
                j = k
                s = 0
        return

    def _calc_motif_undir(self, positions, degrees, k, result, test):
        fn = self.f
        sqn = self.sqsum_r
        for i in range(k):
            if self.r[positions[i]] < degrees[i]:
                return

            fn -= self.r[positions[i]]
            sqn = sqn - (self.r[positions[i]] * self.r[positions[i]] - self.r[positions[i]])
        olda = self.sqsum_r / (2 * self.f)

        newa = sqn / (2 * fn)
        oldb = 0
        be = 0
        for i in range(k):
            for j in range(i + 1, k):
                be += self.r[positions[i]] * self.r[positions[i]]
        newb = be / fn
        changelog = 0
        for i in range(k):
            changelog += self.logfac_r[positions[i]] - math.lgamma(self.r[positions[i]] - degrees[i] + 1)

        x = ((fn * math.log(fn / EXP)) - self.f * math.log(self.f / EXP)) / 2 + olda * olda + olda - newa * \
            newa - newa - newb + changelog

        result.append(x)
        return

    def _calc_motif_dir(self, positions, degrees_r, degrees_c, k, result):
        delta_f = 0
        sqn_r = self.sqsum_r
        sqn_c = self.sqsum_c
        for i in range(k):
            if self.r[positions[i]] < degrees_r[i] or self.c[positions[i]] < degrees_c[i]:
                return
            delta_f += degrees_r[i]
            rold = self.r[positions[i]]
            rnew = self.r[positions[i]] - degrees_r[i]
            sqn_r = sqn_r - (rold * rold - rold) + (rnew * rnew - rnew)
            cold = self.c[positions[i]]
            cnew = self.c[positions[i]] - degrees_c[i]
            sqn_c = sqn_c - (cold * cold - cold) + (cnew * cnew - cnew)

        olda = 0.5 * (self.sqsum_r / self.f) * (self.sqsum_c / self.f)
        newa = 0.5 * (sqn_r / (self.f - delta_f)) * (sqn_c / (self.f - delta_f))
        oldb = self.b
        be = 0
        for i in range(k):
            for j in range(k):
                be += (self.r[positions[i]] - degrees_r[i]) * (self.c[positions[j]] - degrees_c[i])
                if i == j:
                    be -= self.r[positions[i]] * self.c[positions[j]]
        newb = (oldb * self.f + be) / (self.f - delta_f)

        changelog = 0
        for i in range(k):
            changelog += self.logfac_r[positions[i]] - self.logfac_cache[self.r[positions[i]] - degrees_r[i]] + \
                         self.logfac_c[positions[i]] - self.logfac_cache[self.c[positions[i]] - degrees_c[i]]
        x = self.logfac_f[delta_f] - self.logfac_f[0] + olda + oldb - newa - newb + changelog
        result.append(x)
        return


def new_egde(u, v):
    return u << 32 | v


def edge_code(u, v):
    if u < v:
        return u << 32 | v
    else:
        return v << 32 | u


def edge_get_u(e):
    return e >> 32


def edge_get_v(e):
    return e & 4294967295


def reverse(et):
    return (et >> 1) | ((et << 1) & 2)


def DEL(g, row, col):
    g &= ~(1 << (63 - (row * 8 + col)))
    return g


def SET(g, row, col):
    g |= (1 << (63 - (row * 8 + col)))
    return g


def adj(g, k):
    for i in range(k):
        for j in range(k):
            shift = 63 - i * 8 - j
            print(((g>>shift) &1))
        print()
    print()
    return
