import igraph as ig
import numpy as np
import pandas as pd

class tm_ReturnValue:
    def __init__(self, y0, y1, y2):
        self.graph = y0
        self.dataframe = y1
        self.bipartite = y2
        
def k2h_diff(x, d=1):
    return np.array([x[i]-x[i-d] for i in np.arange(d,len(x))])

def K2h(df=None, extras=[], dt=0, dT=0, h=2, bg=None, verbose=False, return_df=False, return_bg=True):
    """
    Finds temporal motifs given dataframe 'df' with required columns ['event_id','actor_id','time','root_id'].
    
    ----------

    df: pandas DataFrame with the required columns ['event_id','actor_id','time','root_id']
        not required if 'bg' is supplied
    extras: list of extra columns to extract from 'df' and pass as attributes

    dt: reaction time for the motifs, in seconds
    dT: repetition time for the motifs, in seconds
    h: the 'h' in K2h, the number of distinct root events in the motifs (>=2)

    bg: bipartite graph obtained from a previous call of the 'K2h' function with the same dataframe 'df';
        if supplied, 'df' is not required

    verbose: if set, return the tuples of root id's for each temporal motif
    return_df: if set, also output the reults as a pandas DataFrame
    return_bg: if set, output the bipartite graph that can be re-used for faster subsequent runs of 'K2h' 
    
    Returns
    -------
        
    self.graph: igraph Graph with K2h motifs (nodes=actors, edges=temporal motifs)
    self.df: pandas DataFrame with K2h motifs (optional, default=False)
    self.bg: igraph multi-bipartite Graph; can be used for faster subsequent runs of 'K2h' (optional, default=True)

    Example
    -------
    Given dataframe D with required columns ['event_id','actor_id','time','root_id']
    
    >> motifs = K2h(df=D, return_df=True, verbose=True, dt=300, dT=1800, h=2)
    >> motifs.df.head() ## output in DataFrame
    >> motifs.graph ## output in igraph Graph
    >> ## re-run using bipartite graph from previous run
    >> new_motifs = K2h(return_df=True, verbose=True, dt=100, dT=1000, h=2, bg=motifs.bipartite)
    
    """
    if bg == None:
        bg = ig.Graph(directed=True)
        src = list(set(df['actor_id']))
        dst = list(set(df['root_id']))
        bg.add_vertices(src)
        bg.add_vertices(dst)
        src = list(df['actor_id'])
        dst = list(df['root_id'])
        edg = [(src[i],dst[i]) for i in range(len(src))]
        bg.add_edges(edg)
        bg.es['time'] = list(df['time'])
        bg.es['is_root'] = list(df['event_id'] == df['root_id'])
        if len(extras)>0:
            for x in extras:
                bg.es[x] = list(df[x])
                
    ## K2h part
    if h<2:
        print('need to have h >= 2')
        return -1
    L = []
    bg.vs['in'] = bg.degree(mode='in')
    ## build all edges from event to root event given dt
    min_time = min(bg.es['time'])
    for v in bg.vs:
        if v['in']>1: ## for each root event tree
            inc = v.incident(mode='in') ## all actors with an event under this root event
            t = [e['time'] for e in inc if e['is_root']] ## time of root event
            if len(t)>0: ## only consider events with a root
                t_min = t[0] ## only 1 root is assumed
                sub = [bg.vs[e.tuple[0]]['name'] for e in inc if e['is_root']][0] ## root actor
                E = list(set([(bg.vs[e.tuple[0]]['name'],sub,t_min-min_time,v['name']) 
                            for e in inc if e['time']-t_min <= dt 
                            and not e['is_root'] ## redundant but slightly faster
                            and sub != bg.vs[e.tuple[0]]['name']])) ## subsequent actors
                l = len(E)
                if l>0:
                    L.extend(E) ## add to list
    g = ig.Graph.TupleList(L, directed=True, edge_attrs=['time','root']) ## build new actor-actor digraph  
    g.es['weight'] = 1
    g.es['time'] = [str(e['time'])+' ' for e in g.es]
    g.es['root'] = [str(e['root'])+' ' for e in g.es]

    ## simplify this graph: add weights, concatenate times and roots
    g = g.simplify(combine_edges=dict(weight=sum, time='concat', root='concat'))
    g.es.select(weight_lt=h).delete() # drop edges with weight < h
    g.vs['deg'] = g.degree() 
    g.vs.select(deg_eq=0).delete() # drop v with degree 0
    
    ## get groups of size h happening within dT
    E = [] ## list of directed edges
    for e in g.es:
        k = 0
        tm = [int(i) for i in e['time'][:-1].split(' ')]
        od = np.argsort(tm) ## times are not sorted
        tm = [tm[i] for i in od]
        th = [i for i in e['root'][:-1].split(' ')]
        th = [th[i] for i in od]
        df = k2h_diff(tm, d=h-1) ## difference between first and last of the h
        trd = []
        
        for j in range(len(df)):
            if df[j]<=dT:
                k+=1
                if verbose:
                    x = []
                    for i in range(h):
                        x.append(th[j+i])
                    trd.append(tuple(x)) ## store tuple of h root events
        if k>0:
            ## store: actor, root actor, count, list of root_id tuples
            if verbose:
                E.append((g.vs[e.tuple[0]]['name'], g.vs[e.tuple[1]]['name'], k, trd)) 
            else:
                E.append((g.vs[e.tuple[0]]['name'], g.vs[e.tuple[1]]['name'], k))
    if verbose:
        G = ig.Graph.TupleList(E, directed=True, edge_attrs=['weight','events'])
    else:
        G = ig.Graph.TupleList(E, directed=True, edge_attrs=['weight'])
    G.vs['in'] = G.degree(mode='IN')
    G.vs['out'] = G.degree(mode='OUT')
    
    ## return DataFrame?
    if return_df:
        G.vs['cc'] = G.clusters(mode='weak').membership
        L = []
        if verbose:
            for e in G.es:
                L.append([G.vs[e.tuple[0]]['cc'],G.vs[e.tuple[0]]['name'],
                          G.vs[e.tuple[1]]['name'],e['weight'],e['events']])
            D = pd.DataFrame(L, columns=['component','root','actor','count','list'])
        else:
            for e in G.es:
                L.append([G.vs[e.tuple[0]]['cc'],G.vs[e.tuple[0]]['name'],
                          G.vs[e.tuple[1]]['name'],e['weight']])
            D = pd.DataFrame(L, columns=['component','root','actor','count'])        
        D.sort_values(by='component', inplace=True)
        D.reset_index(inplace=True, drop=True)
    else:
        D = None
    return tm_ReturnValue(G, D, bg)

## prune G w.r.t. edge weight and color nodes w.r.t. role:
def prune_and_color(g, min_weight=1, min_size=2):    
    G = g.copy()
    ## prune w.r.t. weight 
    e = list(np.where(np.array(G.es['weight'])<min_weight)[0])
    G.delete_edges(e)
    v = list(np.where(np.array(G.degree())==0)[0])
    G.delete_vertices(v)
    ## prune w.r.t. component size
    cls = G.clusters(mode='weak').sizes()
    clm = G.clusters(mode='weak').membership
    G.vs['cls'] = [cls[i] for i in clm]
    G.delete_vertices([v for v in G.vs if v['cls']<=min_size])
    del(G.vs['cls'])
    ## edge width and node colour
    G.es['width'] = [1+np.log2(i) for i in G.es['weight']]
    for v in G.vs:
        if G.degree(v,'in')==0:
            v['color'] = 'blue'
        else:
            if G.degree(v,'out')==0:
                v['color'] = 'red'
            else:
                v['color'] = 'gold'
    return G
