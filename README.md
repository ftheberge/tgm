# tgm - Temporal Graph Motifs

Python code to look for temporal motifs, namely temporal K_{2,h} motifs, in a graph or network of **events**.

We assume that events are grouped under some **root event** which happens first. Some examples are:

* re-tweets given some original tweet (the root event)
* comments following a Reddit submissiion (the root event)

## Data Format

The data consists of temporal 4- or 5-tuples:
* **(event_id, actor_id, time, root_id)**
* **(event_id, actor_id, time, root_id, parent_id)**

Where:
* **event_id**: is a unique id for each event
* **actor_id**: is the id of the actor/user generating the event
* **time**: in UTC format
* **root_id**: is the id of the root event that this event_id is part of
* **parent_id**: is the id of the parent of this event when events with same root_id have a tree structure

Assumption:
* event_id == root_id for the **root event**; in this case, parent_id is meaningless
* all events with same root_id can not happen before the root event itself (where event_id == root_id)
* if "parent_id" are used, the parent event can not happen before one of its child event

## Temporal K_{2,h} motifs

Such motifs are parametrized by 3 values:
* dt: the reaction time, in seconds
* dT: the repetition time, in seconds
* h: the number of distinct root events forming the motif

A temporal K_{2,h} motifs given (dt, dT) occurs when:
* an actor A submits a root event and a different actor B sumbits an event under that root event (i.e. same root_id) within dt seconds (**reaction** time)
* the above scenario happens h times, for h distinct root events, withith dT seconds (**repetition** time)

## Finding K_{2,h} motifs

The function **K2h** is called to find temporal motifs. When called first, a DataFrame with the required columns needs to be passed (see details below). If the function is called several times for the same dataset (for example, trying different dt and dT parameters), then some speed-up can be achieved as we illustrate below.

### Parameters:

* df: pandas DataFrame with the required columns:
  * **event_id**: unique id for each event
  * **author_id**: actor/user id
  * **time**: in seconds, typically UTC
  * **root_id**: id of the root event that this event_id is part of
* extras: list of extra columns to extract from 'df' and keep as attributes
* dt: reaction time for the motifs, in seconds
* dT: repetition time for the motifs, in seconds
* h: the 'h' in K2h; the number of different root_id's in the motifs (>=2)
* bg: bipartite graph obtained from a previous call of the 'K2h' function with the same dataframe 'df'; if supplied, 'df' is not required; this will speed-up the computation
* verbose: if set, return tuples of root_id's for each motif (default=False)
* return_df: if set, also output the results as a pandas DataFrame (default=False)
* return_bg: if set, output the bipartite graph that can be re-used for faster subsequent runs of 'K2h' (default=True)

Result are pairs of actors with the number of temporal motifs they share and if 'verbose' is set, the list of tuples of root_id's they share. 

### Output:    
* self.graph: igraph Graph with the K2h motifs (nodes=actors, edges=K_{2,h} motifs)
* self.df: pandas DataFrame with K2h motifs (optional, default=False)
* self.bg: igraph multi-bipartite Graph; can be used for faster subsequent runs of 'K2h' (optional, default=True)

For the DataFrame output, we also group the actors into weakly connected components, returning the component number for each one.

For the Graph output, the edges are **directed** from the non-root actor to the root actor. Edge weights indicate the number of motifs. There are no loops.
