#!/usr/bin/python

pdhelp="""Process data.

   Run unit tests:
     ./pd.py -t
     
   Parse KokkosKernels_Test_{routine}* files and do some simple plots:
     ./pd.py --kkt-parse --kkt-plot-vs-nthreads-loglog -f KokkosKernels_Test_Gemm.txt
     ./pd.py --kkt-parse --kkt-plot-vs-nthreads-loglog -f KokkosKernels_Test_LU.txt
     ./pd.py --kkt-parse --kkt-plot-vs-nthreads-loglog -f KokkosKernels_Test_Trsm.txt
     
   Save the figure to pdf instead of showing it:
     ./pd.py --kkt-parse --kkt-plot-vs-nthreads-loglog -f KokkosKernels_Test_Gemm.txt -s foo
     
   Linear-linear:
     ./pd.py --kkt-parse --kkt-plot-vs-nthreads-linlin -f KokkosKernels_Test_Trsm.txt
   Linear-log:
     ./pd.py --kkt-parse --kkt-plot-vs-nthreads-linlog -f KokkosKernels_Test_Trsm.txt

   Show data per thread rather than in total: Add --per-thread.

   Right now, this is probably the best format:
     ./pd.py -t --kkt-parse --kkt-plot-vs-nthreads-linlog -s foo -f KokkosKernels_Test_Gemm.txt --per-thread

   Plots for workset size:
     ./pd.py --kkt-parse --kkt-plot-workset -f KokkosKernels_Test_Gemm.workset.txt -s foo
   To show speedup w.r.t. OpenMP MKL, add --speedup:
     ./pd.py --kkt-parse --kkt-plot-workset --speedup -f KokkosKernels_Test_Gemm.workset.txt -s foo

   To make wide rather than tall figures, add
     --wide
"""

import os, sys, re, optparse
import matplotlib.pyplot as pl
import numpy as np

#> Utils.

def expect(name, a, b):
    'Unit test wrapper.'
    if a != b:
        print '{0} failed: value = {1}; expected = {2}'.format(name, a, b)

def dispfig(fn_prefix):
    pl.tight_layout()
    if len(fn_prefix) == 0:
        pl.show()
    else:
        pl.savefig(fn_prefix + '.pdf', format='pdf', bbox_inches='tight')

def readall(fn):
    'Shorthand for reading in all the text in a file.'
    try:
        with open(fn, 'r') as f:
            text = f.read()
    except:
        text = ''
    return text

def sscanf(string, format):
    """Kind of like sscanf. Format is like this: 'i,if,f,s', where if is for
    int(float(.))."""
    str2type = {'i': int, 'if': lambda x: int(float(x)), 'f': float, 's': str}
    ts = [str2type[s] for s in format.split(',')]
    v = [t(s) for t,s in zip(ts, string.split())]
    return v

def get_first_word(ln):
    'Get first word from a line.'
    s = ln.strip().split()
    if len(s) == 0: return ''
    return s[0]

def get_last_word(ln):
    s = ln.strip().split()
    if len(s) == 0: return ''
    return s[-1]

def strleq(s, ref):
    if len(s) < len(ref): return False;
    return s[0:len(ref)] == ref

def flop_lu(n):
    "# of + and * for LU of nxn matrix."
    return 2.0/3*(n-1)**3 + (n-1)**2 + 1.0/3*(n-1) + n**2/2.0 + n/2.0

def flop_bothtrsm(n, k):
    """# of + and * for both lower and upper trsms, one with unit diag, of nxn
    matrix with k vectors."""
    return (2*n**2 - n)*k

def flop_gemm(n, k):
    """# of + and * for matmat of nxn matrix with nxk matrix, with accumulation
    into the output."""
    return 2*n**2*k

def flop_tridiag_factor(n, b):
    """#FLOP for block tridiag factorization. Block size is b. There are n block
    rows."""
    return n*flop_lu(b) + (n-1)*flop_bothtrsm(b, b) + (n-1)*flop_gemm(b, b)

def flop_tridiag_solve(n, b, k):
    """#FLOP for block tridiag solve. Block size is b. There are n block
    rows. #RHS is k."""
    return n*flop_bothtrsm(b, 1) + 2*(n-1)*flop_gemm(b, 1)

def rename(name):
    "Possibly change the name of something."
    names = {'SIMD': 'PP-Open-Src',
             'Serial SIMD': 'PP-Open-Src',
             'Vector': 'PP-Open-Src',
             'MKL DGEMM': 'MKL OpenMP',
             'MKL Cmpct': 'Vendor-Compact',
             'MKL LU': 'MKL OpenMP',
             'MKL TRSM': 'MKL OpenMP',
             'MKL Cmpt': 'Vendor-Compact',
             'MKLCmpt::': 'Vendor-Inner-Compact',
             'MKLBatch::': 'Vendor-Outer-Compact'}
    if names.has_key(name): return names[name]
    return name

def plot_styles():
    return {'markersize': 8, 'markeredgecolor': 'k', 'markerfacecolor': 'k'}

def patternmap():
    return {'MKL Batch': 'bo-', 'Vendor-Compact': 'r*-', 'MKL OpenMP': 'gs-', 'PP-Open-Src': 'kd-', 'libxsmm': 'cv-',
            'Vendor-Outer-Compact': 'bo-', 'Vendor-Inner-Compact': 'r*-', 'Native': 'gs-'}

#> BLAS/LAPACK analysis.

def kkt_is_data_oneliner(ln):
    'Is this line a data line?'
    if len(ln) < 96:
        return False
    fw = get_first_word(ln)
    return fw == 'MKL' or fw == 'Serial' or fw == 'SIMD' or fw == 'libxsmm' ## or fw == 'Plain' or fw == 'Team' or

def kkt_is_exespace_line(ln):
    'Is this a line giving the execution space and its parameters?'
    if len(ln) < 10:
        return False
    fw = get_first_word(ln)
    return len(fw) >= 9 and fw[0:9] == 'ExecSpace'

def kkt_parse_oneliner_re(ln):
    'Parse a oneliner to a tuple of strings.'
    hits = re.findall('(?P<name>.*) BlkSize = (?P<blksz>[^a-zA-Z]*) .*time = (?P<time>.*) ' + \
                      'avg flop/s = (?P<avgflops>.*) max flop/s = (?P<maxflops>[^ds]*) ', ln + ' ')
    if len(hits) == 0:
        return ()
    return hits[0]

def kkt_parse_oneliner(ln):
    'Parse one line of data -> (alg name, block size, time, avg flop/s, max flop/s).'
    v = kkt_parse_oneliner_re(ln)
    if len(v) == 0:
        return ()
    if len(v) != 5:
        print 'Parse error; ln = ' + ln
        return ()
    return (v[0].strip(), int(v[1]), float(v[2]), float(v[3]), float(v[4]))

def kkt_parse_exespace_re(ln):
    'Parse an ExecSpace line.'
    hits = re.findall('topology\[(?P<one>.*) x (?P<two>.*) x (?P<three>.*) \]', ln + ' ')
    if len(hits) == 0:
        return ()
    return hits[0]

def kkt_parse_exespace(ln):
    'Parse an ExecSpace line -> (#nodes, #cores, #threads/core).'
    v = kkt_parse_exespace_re(ln)
    if len(v) == 0:
        return ()
    if len(v) != 3:
        print 'Parse error; ln = ' + ln
        return ()
    return tuple([int(i) for i in v])

def kkt_is_workset_line(ln): return strleq(ln, ' N = ')
def kkt_parse_workset(ln): return int(ln[5:])

def kkt_parse(text, want_workset):
    "Parse output from K.K.'s performance test driver."
    lns = text.split('\n')
    d = {}
    for ln in lns:
        if kkt_is_exespace_line(ln):
            exespace = kkt_parse_exespace(ln)
            if not d.has_key(exespace):
                d[exespace] = {}
        elif want_workset and kkt_is_workset_line(ln):
            workset = kkt_parse_workset(ln)
        elif kkt_is_data_oneliner(ln):
            v = kkt_parse_oneliner(ln)
            alg = v[0]
            blksz = v[1]
            if not d[exespace].has_key(alg): d[exespace][alg] = {}
            if want_workset:
                if not d[exespace][alg].has_key(blksz): d[exespace][alg][blksz] = {}
                d[exespace][alg][blksz][workset] = v[2:]
            else:
                d[exespace][alg][blksz] = v[2:]
    return d

def kkt_get_nthreads(d):
    """Return sorted list of #threads used in the data. Return also the
    corresponding exespace tuples."""
    s = set()
    for es in d.keys():
        s.add((es[1]*es[2], es))
    s = list(s)
    s.sort()
    return [i[0] for i in s], [i[1] for i in s]

def kkt_get_algorithms(d):
    'Return list of algorithm names.'
    a = d.values()[0].keys()
    a.sort()
    return a

def kkt_get_blocksizes(d):
    'Return list of block sizes.'
    a = d.items()[0][1].items()[0][1].keys()
    a.sort()
    return a

def kkt_plot_vs_nthreads(d, fn_prefix, xlinear=False, ylinear=False, perthread=False, wide=False):
    nthreads, exespaces = kkt_get_nthreads(d)
    algs = kkt_get_algorithms(d)
    blkszs = kkt_get_blocksizes(d)

    corecnt = max(nthreads) == 68
    fs = 12 if wide else 15
    ival = 1
    scale = 1e-9
    x = [1, 2, 4, 8, 16, 34, 68, 136, 272]
    if xlinear:
        x.remove(2)
        xfn = lambda x: x
    else:
        xfn = np.log2
    patmap = patternmap()
    pl.figure(num=1, figsize=(16,4) if wide else (7,16))
    for ibs, bs in enumerate(blkszs):
        if wide: pl.subplot(1, len(blkszs), ibs+1)
        else: pl.subplot(len(blkszs), 1, ibs+1)
        ymin = 10
        ymax = 0
        for ialg, alg in enumerate(algs):
            y = []
            pat = patmap[rename(alg)] if patmap.has_key(rename(alg)) else ''
            for i, nthr in enumerate(nthreads):
                f = scale
                if perthread: 
                    f /= nthr
                y.append(d[exespaces[i]][alg][bs][ival] * f)

            ymin = min(min(y), ymin)
            ymax = max(max(y), ymax)

            if ylinear:
                pl.plot(xfn(nthreads), y, pat, linewidth=2.0, label=rename(alg), **plot_styles())
            else:
                pl.semilogy(xfn(nthreads), y, pat, linewidth=2.0, label=rename(alg), **plot_styles())
        pl.xticks(xfn(x), [str(i) for i in x], fontsize=fs)
        pl.yticks(fontsize=fs)
        ymin *= 0.96
        if ylinear: ymin = 0
        ymax *= 1.04
        xmin = 1;
        xmax = nthreads[-1]
        if xlinear:
            xmin -= 1
            xmax += 1
        else:
            xmin = -0.1
            xmax = np.log2(xmax) + 0.1
        pl.axis([xmin, xmax, ymin, ymax])
        pl.grid(True)

        perthread_str = ''
        if perthread:
            if corecnt: perthread_str = ' / Core'
            else: perthread_str = ' / Thread'
        if ibs == 0 or not wide:
            pl.ylabel('GFLOPS{0}'.format(perthread_str), fontsize=fs)
        pl.title('Block Size {0}'.format(bs))
        if ibs == 0:
            pl.legend(loc='upper right' if perthread else 'upper left', fontsize=fs-1)

        if ibs == 3 or wide:
            pl.xlabel('# Cores' if corecnt else '# Threads', fontsize=fs)

    dispfig(fn_prefix)

def plot_lol(lol, xlabel=None, ylabel=None, title=None, fontsize=20, xvals=None, yvals=None, scale=1.0):
    """lol is a list of lists indexed as lol[ix][yy]. The y-direction lists should
    all be the same size."""
    ax = pl.gca()
    nx = len(lol)
    ny = len(lol[0])
    xsz = 1.0/nx
    ysz = 1.0/ny
    for i in range(ny):
        for j in range(nx):
            xll = j/float(nx) - xsz/2
            yll = i/float(ny) - ysz/2
            r = pl.Rectangle([xll, yll], xsz, ysz, facecolor='black', edgecolor='white',
                             alpha=min(0.6, float(lol[j][i])*scale))
            ax.add_patch(r)
            pl.text(xll + xsz/2, yll + ysz/2, '{0}'.format(lol[j][i]),
                    horizontalalignment='center', verticalalignment='center', fontsize=fontsize)
    pl.axis([-xsz/2, (nx-0.5)*xsz, -ysz/2, (ny-0.5)*ysz])
    xticks = [float(x)/nx for x in range(nx)]
    yticks = [float(y)/ny for y in range(ny)]
    pl.xticks(xticks, xvals, fontsize=fontsize)
    pl.yticks(yticks, yvals, fontsize=fontsize)
    if xlabel: pl.xlabel(xlabel, fontsize=fontsize)
    if ylabel: pl.ylabel(ylabel, fontsize=fontsize)
    if title: pl.title(title, fontsize=fontsize)

def test_plot_lol(fn_prefix):
    lol = []
    for j in range(5):
        lol.append([])
        for i in range(4):
            lol[j].append((j,i))
    plot_lol(lol, xlabel='xlbl', ylabel='ylbl', title='ttl')
    dispfig(fn_prefix)

def kkt_plot_workset(d, fn_prefix, speedup=False, wide=False):
    assert(len(d.keys()) == 1)
    d0 = d.values()[0]
    
    algs = d0.keys()
    algs.sort()
    blkszs = d0.values()[0].keys()
    blkszs.sort()
    wrksets = d0.values()[0].values()[0].keys()
    wrksets.sort()
    ival = 1

    if speedup:
        colorscale = 1.0/40
        for n in ['MKL LU', 'MKL DGEMM', 'MKL TRSM']:
            try:
                algs.remove(n)
                against = n
            except:
                pass
    else:
        colorscale = 1.0/500
        # Put OpenMP at the end.
        omp = 'MKL OpenMP'
        try:
            idx = [rename(a) for a in algs].index(omp)
            print idx
            algs = algs[:idx] + algs[idx+1:] + [algs[idx]]
            print algs
        except: pass

    # Don't show KK in these plots.
    for a in ['Serial SIMD', 'SIMD']:
        try: algs.remove(a)
        except: pass

    pl.figure(num=1, figsize=((4 if speedup else 4)*len(algs), 3) if wide else (5, 3*len(algs)))
    scale = 1e-9
    aspect = float(len(blkszs))/len(wrksets)
    for ia, alg in enumerate(algs):
        if wide: pl.subplot(1, len(algs), ia+1, aspect=aspect)
        else: pl.subplot(len(algs), 1, ia+1, aspect=aspect)
        lol = []
        for ws in wrksets:
            y = []
            for bs in blkszs:
                if speedup: scale = 1.0 / d0[against][bs][ws][ival]
                y.append('{0:1.1f}'.format(d0[alg][bs][ws][ival]*scale))
            lol.append(y)
        plot_lol(lol, fontsize=11,
                 xlabel=r'log$_2$ Workset Size' if wide or ia == len(algs)-1 else None,
                 ylabel='Block Size' if not wide or ia == 0 else None,
                 title=rename(alg),
                 xvals=[int(round(np.log2(w))) for w in wrksets], yvals=blkszs, scale=colorscale)

    dispfig(fn_prefix)

#> Driver.

def run_tests():
    'Unit tests for the above.'
    expect('scanf', sscanf('  hi 3.14e-1  7', 's,f,i'), ['hi', 0.314, 7])
    expect('get_first_word', get_first_word('   MKL jumped  over = '), 'MKL')
    expect('kkt_is_data_oneliner', kkt_is_data_oneliner('moo'), False)
    testln1 = '  MKL Cmpt BlkSize =   9 time = 4.460000e-04 avg flop/s = 4.237677e+09 max flop/s = 1.664115e+10 diff to ref = 0.000000e+00'
    testln2 = '  MKL Cmpt BlkSize =   9 time = 4.460000e-04 avg flop/s = 4.237677e+09 max flop/s = 1.664115e+10'
    testln3 = 'ExecSpace::  Kokkos::OpenMP KOKKOS_ENABLE_OPENMP hwloc[1x68x4] hwloc_binding_enabled thread_pool_topology[ 1 x 68 x 4 ]'
    expect('kkt_is_data_oneliner', kkt_is_data_oneliner(testln1), True)
    expect('kkt_parse_oneliner_re', kkt_parse_oneliner_re(testln2), ('  MKL Cmpt', '  9', '4.460000e-04', '4.237677e+09', '1.664115e+10'))
    expect('kkt_parse_oneliner_re', kkt_parse_oneliner_re(testln1), ('  MKL Cmpt', '  9', '4.460000e-04', '4.237677e+09', '1.664115e+10'))
    expect('kkt_parse_exespace', kkt_parse_exespace(testln3), (1, 68, 4))
    expect('flop_lu', flop_lu(17), 3145.0)
    expect('flop_bothtrsm', flop_bothtrsm(5, 4), 180)
    expect('flop_gemm', flop_gemm(5, 4), 200)

def get_optparser():
    'Command-line parser.'
    p = optparse.OptionParser(usage=pdhelp)
    p.add_option('-t', '--test', dest='test', action='store_true', default=False, help='Run unit tests.')
    p.add_option('-f', '--filename', dest='filename', default='', help='Filename to parse.')
    
    p.add_option('--kkt-parse', dest='kkt_parse', action='store_true', default=False, help='Parse KokkosKernels_Test_{routine}* files.')
    p.add_option('--kkt-plot-vs-nthreads-loglog', dest='kkt_plot_vs_nthreads_loglog', action='store_true', default=False)
    p.add_option('--kkt-plot-vs-nthreads-linlin', dest='kkt_plot_vs_nthreads_linlin', action='store_true', default=False)
    p.add_option('--kkt-plot-vs-nthreads-linlog', dest='kkt_plot_vs_nthreads_linlog', action='store_true', default=False)
    p.add_option('--kkt-plot-workset', dest='kkt_plot_workset', action='store_true', default=False)

    p.add_option('--per-thread', dest='perthread', action='store_true', default=False)
    p.add_option('--gflops', dest='gflops', action='store_true', default=False)
    p.add_option('--speedup', dest='speedup', action='store_true', default=False)
    p.add_option('--wide', dest='wide', action='store_true', default=False)
    
    p.add_option('-s', '--save', dest='save_prefix', default='')
    return p

if __name__ == '__main__':
    p = get_optparser()
    (opts, args) = p.parse_args()
    if opts.test: run_tests()

    if opts.kkt_parse:
        text = readall(opts.filename)
        if len(text) == 0:
            print "Empty file '" + opts.filename + "'; exiting."
            sys.exit(-1)        

    if opts.kkt_parse:
        d = kkt_parse(text, want_workset=opts.kkt_plot_workset)
        if opts.kkt_plot_vs_nthreads_loglog:
            kkt_plot_vs_nthreads(d, opts.save_prefix, xlinear=False, ylinear=False,
                                 perthread=opts.perthread, wide=opts.wide)
        elif opts.kkt_plot_vs_nthreads_linlin:
            kkt_plot_vs_nthreads(d, opts.save_prefix, xlinear=True, ylinear=True,
                                 perthread=opts.perthread, wide=opts.wide)
        elif opts.kkt_plot_vs_nthreads_linlog:
            kkt_plot_vs_nthreads(d, opts.save_prefix, xlinear=False, ylinear=True,
                                 perthread=opts.perthread, wide=opts.wide)
        elif opts.kkt_plot_workset:
            kkt_plot_workset(d, opts.save_prefix, speedup=opts.speedup, wide=opts.wide)
