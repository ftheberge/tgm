# tgm - Temporal Graph Motifs

Python code to look for temporal motifs, namely temporal K_{2,h} motifs, in a network. 

## Data Format

The data consists of temporal 4- or 5-tuples:
* **(event_id, author_id, time, root_id)**
* **(event_id, author_id, time, root_id, parent_id)**

Where:
* event_id: is a unique id for each event
* author_id: is the id of the author/user generating the event
* time: in UTC format
* root_id: is the id of the root event that this event_id is part of
* parent_id: is the id of the parent of this event when events with same root_id have a tree structure

Assumption:
* event_id == root_id for the root event; in this case, parent_id is meaningless
* all events with same root_id can not happen before the root event itself (where event_id == root_id)
* if "parent_id" are used, the parent event can not h

## Finding K_{2,h} motifs

The function **K2h** is called to find temporal motifs. When called first, a DataFrame with the required columns needs to be passed. If the function is called several times for the same dataset (for example, trying different dt and dT parameters), then some seepd-up can be achieved as we illustrate below.

### Parameters:

* df: pandas DataFrame with the required columns:
 * **event_id**: unique id for each event
 * **author_id**: a.k.a. userid
 * **time**: in seconds, typically UTC
 * **root_id**: id of the root event that this event_id is part of
* extras: list of extra columns to extract from 'df' and keep as attributes
* dt: reaction time for the motifs, in seconds
* dT: repetition time for the motifs, in seconds
* h: the 'h' in K2h, the number of different root_id's in the motifs (>=2)
* bg: bipartite graph obtained from a previous call of the 'K2h' function with the same dataframe 'df'; if supplied, 'df' is not required; this will speed-up the computation
* verbose: if set, return tuples of root_id's for each motif (default=False)
* return_df: if set, also output the results as a pandas DataFrame (default=False)
* return_bg: if set, output the bipartite graph that can be re-used for faster subsequent runs of 'K2h' (default=True)

Result are pairs of authors, the number of temporal motifs they share and if 'verbose' is set, the list of tuples of root_id's they share. 

### Output:    
* self.graph: igraph Graph with K2h motifs
* self.df: pandas DataFrame with K2h motifs (optional, default=False)
* self.bg: igraph multi-bipartite Graph; can be used for faster subsequent runs of 'K2h' (optional, default=True)

For the DataFrame output, we also group the authors into weakly connected components, returning the component number for each one.


