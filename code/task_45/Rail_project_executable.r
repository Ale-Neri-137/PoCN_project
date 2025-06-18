
library(sf)
library(dplyr)
library(igraph)
library(visNetwork)
library(geosphere)
library(shiny)
library(htmlwidgets)
library(tidyr)     # unnest_longer
library(purrr)     # map2_dfr
library(lwgeom)
r_world <- 6378137






# funzione che data una edge list, restituisce 
# le coordinate della node list in formato dataframe (id, x, y)
assign_coord <- function(coord_dictionary, edge_list){
    
    # Estrai l'insieme unico dei nodi dagli archi
    # questa cosa da problemi nel caso in cui ci siano nodi isolati!!!
    node_list <- unique(as.vector(edge_list))

    # Filtra la lista delle coordinate dei vertici in base ai nodi presenti in node_set

    coord_list <- coord_dictionary[coord_dictionary$id %in% node_list,]
    return(coord_list)
}
# dato un grafo ( e gli id delle stazioni )
# restituisce i nodi del grafo 'ridondanti' e i due neighbors
# redundat nodes are degree 2 nodes that aren't stations, 
# for sure we can remove those and connect the 2 neighbors
find_redundancy_matrix <- function(graph, protected_nodes){

    redundant_id <- V(graph)[!(V(graph)$name %in% protected_nodes) & degree(graph) == 2]$name
    # Initialize a list to store the results
    neighbor_a <- c()
    neighbor_b <- c()
    
    # Loop through each redundant node and get its neighbors
    for(id in redundant_id) {
        #neighbors_id <- neighbors(graph, id)
        neighbor_a <- c(neighbor_a, as.double(neighbors(graph, id)[[1]]$name))
        neighbor_b <- c(neighbor_b, as.double(neighbors(graph, id)[[2]]$name))
    }
    
    redundancy_mat <- data.frame(as.double(redundant_id), neighbor_a, neighbor_b)
    return(as.matrix(redundancy_mat))
}

flip_index <- function(idx){
    if(idx == 1){
        return(2)
    }
    else if(idx == 2){
        return(1)
    }
}

# questa funzione si prende la redundancy matrix, un grafo e restituisce il grafo aggiornato 
# semplificato dall'algoritmo 
simplify_graph_one_step <- function(redundancy_matrix,  graph){

    edgelist <- as_edgelist(graph = graph, names = TRUE)
    redundant_nodes <- redundancy_matrix[,1]
    neighbors_mat <- redundancy_matrix[,2:3, drop = FALSE]#
    new_edge_list <- matrix(ncol = 2, nrow = 0) 

    # Trova e salva shortcut edge tra nodi non ridondanti
    for(i in 1:nrow(neighbors_mat)){
        for(j in 1:ncol(neighbors_mat)){

            #controlliamo se elemento è brutto 
            if((neighbors_mat[i,j] %in% redundant_nodes)){
                next
            }
            #se bello 
            else {
                #salva elemento
                new_edge_list <- rbind(new_edge_list, c(neighbors_mat[i,j],0)) 
                
                # vedi suo amico
                friend <- neighbors_mat[i,flip_index(j)]
                # cerca se esiste in ab il twin
                twin_index <- which(row(neighbors_mat) != i & neighbors_mat == friend, arr.ind = TRUE)
                # verifichiamo che friend sia ridondante e che il suo twin esista 
                condition <- friend %in% redundant_nodes && !is.na(twin_index[1])

                while(condition){
                    # aggiorniamo friend 
                    friend <- neighbors_mat[twin_index[1],flip_index(twin_index[2])]
                    i <- twin_index[1]
                    twin_index <- which(row(neighbors_mat) != i & neighbors_mat == friend, arr.ind = TRUE)
                    condition <- !is.na(twin_index[1]) && friend %in% redundant_nodes 
                }
                
                if(!(friend %in% redundant_nodes)){
                    new_edge_list[nrow(new_edge_list),2] <- friend 
                    
                }
                else if (is.na(twin_index[1])){
                    i_prime <- which(redundant_nodes == friend ,  arr.ind = TRUE)
                    row_elements <- neighbors_mat[i_prime,]
                    element <- row_elements[!row_elements %in% redundant_nodes]
                    new_edge_list[nrow(new_edge_list),2] <- (element)
                }
            }
        }
    }
        # Conserva gli edge tra due nodi non ridondanti
    for(i in 1:nrow(edgelist)){
        row_elements <- edgelist[i,]
        if(!any(row_elements %in% redundant_nodes)){
            new_edge_list <- rbind(new_edge_list, row_elements)
        }
    }
    return(unique(new_edge_list))
}

get_protected_nodes <- function(g, stations_id) {
  # 1) Stations are automatically protected
  prot <- stations_id
  
  # 2) Find all 3-cliques and mark their members
  tris <- cliques(g, min = 3, max = 3)
  tri_members <- unique(unlist(lapply(tris, function(tri) V(g)[tri]$name)))
  prot <- c(prot, tri_members)
  
  # 3) Find all nodes j with degree >= 3 that are NOT in any 3-clique
  degs <- degree(g)
  names(degs) <- V(g)$name
  high_deg <- names(degs)[degs >= 3]
  tri_members_deg3 <- setdiff(tri_members, names(degs)[degs >= 4])
  non_tri_high_deg <- setdiff(high_deg, tri_members_deg3)
  
  # 4) For each such j, include *all* of its neighbors
  nbrs <- unique(unlist(neighborhood(g, order = 1, nodes = non_tri_high_deg, mode = "all")) )
  nbrs <- V(g)[nbrs]$name
  
  # 5) Combine and unique
  unique(c(prot, nbrs))
}

simplify_graph <- function(graph){

    protected_nodes <- get_protected_nodes(g,stations_id)
    non_protected <- setdiff(V(graph)$name, protected_nodes)
    redundancy_matrix <- find_redundancy_matrix(graph, protected_nodes)  

    condition <- any(degree(graph, non_protected) == 2)
    graph_history <- list(graph)
    i <- 1

    while(condition){

        simplified_edge_list <- simplify_graph_one_step(redundancy_matrix, graph)

        vert <- as_data_frame(g, what = "vertices")
        used_names <- unique(c(simplified_edge_list[, 1], simplified_edge_list[, 2])
)      
        vert   <- vert %>%          
         filter(name %in% used_names)

        graph <- igraph::simplify(graph_from_data_frame(simplified_edge_list,
                                                 directed = F,vertices = vert))
        graph_history[[i+1]] <- graph

        protected_nodes <- get_protected_nodes(g,stations_id)
        non_protected <- setdiff(V(graph)$name, protected_nodes)
        redundancy_matrix <- find_redundancy_matrix(graph, protected_nodes)        
        
        condition <- any(degree(graph, non_protected)==2) & nrow(redundancy_matrix)!=0
        i <- i+1
    }
    cat('\n grafo semplificato :) !')
    return(graph_history)   
}

find_contractible_Ys <- function(g, stations_id, maxdist = 3) {
  # 1) candidate junctions: degree-3 & not a station
  degs      <- degree(g)
  candidates <- V(g)[degs == 3]$name
  
  triples <- list()
  
  for (j in candidates) {
    # 2) find other degree-3 nodes within distance ≤ maxdist
    others      <- setdiff(names(degs)[degs == 3], j)
    d_to_others <- distances(g, v = j, to = others)[1, ]
    close_nodes <- others[d_to_others <= maxdist]
    if (length(close_nodes) < 2) next
    
    # 3) examine every pair (k1,k2) among close_nodes
    for (pair in combn(close_nodes, 2, simplify = FALSE)) {
      k1 <- pair[1]; k2 <- pair[2]
      
      # 4) check pairwise distance
      if (distances(g, v = k1, to = k2)[1,1] > maxdist) next
      
      # 5) get the three shortest‐path vertex sequences
      paths <- list(
        p1 = shortest_paths(g, from = j,  to = k1, output = "vpath")$vpath[[1]],
        p2 = shortest_paths(g, from = j,  to = k2, output = "vpath")$vpath[[1]],
        p3 = shortest_paths(g, from = k1, to = k2,output = "vpath")$vpath[[1]]
      )
      
      # 6) verify each intermediate node on each path
      ok <- TRUE
      for (p in paths) {
        verts <- as_ids(p)
        if (length(verts) > 2) {
          mids <- verts[-c(1, length(verts))]
          # must not be a station, and must have degree 2
          if (any(mids %in% stations_id) ||
              any(degree(g, mids) != 2)) {
            ok <- FALSE
            break
          }
        }
      }
      if (ok) {
        triples[[length(triples) + 1]] <- c(j, k1, k2)
      }
    }
  }
  
  # dedupe and return
  unique(triples)
}

collapse_contractible_Ys <- function(g, stations_id, maxdist = 3) {
  triples <- find_contractible_Ys(g, stations_id, maxdist)
  g2      <- g

  for (tri in triples) {
    pairs <- combn(tri, 2, simplify = FALSE)

    # 1) Collect mids on the original g
    mids <- unique(unlist(lapply(pairs, function(p) {
      pth <- shortest_paths(
        g, from = p[1], to = p[2], output = "vpath"
      )$vpath[[1]]
      ids <- as_ids(pth)
      if (length(ids) > 2) ids[-c(1, length(ids))] else character(0)
    })))

    # 2) Add the missing clique‐edges to g2
    for (p in pairs) {
      if (!are.connected(g2, p[1], p[2])) {
        g2 <- add_edges(g2, c(p[1], p[2]))
      }
    }

    # 3) Only delete those mids that still exist in g2
    mids_to_rm <- intersect(mids, V(g2)$name)
    if (length(mids_to_rm)) {
      g2 <- delete_vertices(g2, mids_to_rm)
    }
  }

  # 4) final clean‐up
  igraph::simplify(g2)
}

##Now data extraction and graph simplification starts, in one giant loop over countries

countries <- c("AL","AT","BE","BG","CHLI","CY","CZ","DE","DK",
               "EE","ES","FI","FR","GB","GE","GR","HR","HU",
               "IE","IS")

# 2) iterate ---------------------------------------------------------------
for (icc in countries) {

  shp_pts <- file.path(
    "Data_task_45/EGM_2019_SHP_20190312",
    "DATA","Countries", icc, "RailrdC.shp")

  if (!file.exists(shp_pts)) {
    warning("missing file: ", shp_pts)
    next                      # skip to next code
  }

  shp_rail <- file.path(
    "Data_task_45/EGM_2019_SHP_20190312",
    "DATA","Countries", icc, "RailrdL.shp")

  if (!file.exists(shp_rail)) {
    warning("missing file: ", shp_rail)
    next                      # skip to next code
  }


  message("• reading ", icc, " …")

lines        <- st_read(shp_rail, quiet = TRUE) 

# 2) Extract raw coords + feature-index
xy           <- st_coordinates(lines)   # columns: X, Y, L1

# 3) Pull out just the two attrs we want and tag them with L1
attrs_by_L1  <- lines %>%
                  st_drop_geometry() %>%      # drop geometry column
                  mutate(L1 = row_number()) %>% 
                  select(L1, ICC)

# 4) Build coord table and join attributes
coords_rail    <- as_tibble(xy) %>%         # X, Y, L1
                  left_join(attrs_by_L1, by = "L1") %>%
                  group_by(L1) %>%
                  mutate(
                    pos = row_number(),
                    n   = n()
                  ) %>%
                  ungroup()

# 1) Remove the big intermediate objects
rm(lines, xy)  
# 2) Run garbage‐collection
invisible(gc())

coords_rail_pruned <- coords_rail %>%
  filter(pos <= 2     # first two points
        | pos >= n - 1 # last two points
        )

# 2) Read & force to WGS‐84 lon/lat
stations_sf <- st_read(shp_pts, quiet = TRUE) 
#print(st_crs(stations_sf))
#print(unique(stations_sf$TFC))
# 3) Extract raw coords (X=lon, Y=lat)
xy_pts      <- st_coordinates(stations_sf)

# 4) Drop the geometry and grab only the attrs you care about
attrs       <- stations_sf %>%
  st_drop_geometry() %>%
  select(ICC, NAMN1, RStationID,TFC)

# 5) Assemble the stat_coords tibble
coords_stations <- tibble(
  lon          = xy_pts[, "X"],
  lat          = xy_pts[, "Y"]
) %>%
  bind_cols(attrs)

# 1) Remove the big intermediate objects
rm(stations_sf, xy_pts)  
# 2) Run garbage‐collection
invisible(gc())

coords_stations <- coords_stations %>%
  distinct() 

coords_rail_enriched <- coords_rail_pruned %>% #you can use the pruned or not pruned here
  # join on exact coordinate match
  left_join(
    coords_stations,
    by = c("X" = "lon", "Y" = "lat"),
    relationship = "many-to-many"
  )
coords_rail_enriched <- coords_rail_enriched %>%
  # drop the one we don’t want:
  select(-ICC.y) %>%        
  # rename the survivor to a single “ICC”:
  rename(ICC = ICC.x)       

  nodes <- coords_rail_enriched %>%
  group_by(X, Y) %>%
  summarise(
    ICC    = list(sort(unique(ICC))),   # a list‐column of codes, e.g. c("ES","PT")
    NAMN1  = first(na.omit(NAMN1)) %||% NA_character_,
    .groups = "drop"
  ) %>%
  mutate(node_id = row_number())

# 2) Build the edge list ---------------------------------------------------

edges <- coords_rail_enriched %>%
  arrange(L1, pos) %>%         # ensure each L1 is in order
  group_by(L1) %>%
  mutate(
    X_to = lead(X),            # next‐point coords
    Y_to = lead(Y)
  ) %>%
  filter(!is.na(X_to)) %>%     # drop last point of each line (no lead)
  ungroup() %>%
  # now map coords → node_id
  left_join(nodes, by = c("X" = "X", "Y" = "Y")) %>%
  rename(from = node_id) %>%
  left_join(nodes, by = c("X_to" = "X", "Y_to" = "Y")) %>%
  rename(to   = node_id) %>%
  select(from, to)

# 3) Create the igraph object ----------------------------------------------

g <- graph_from_data_frame(
  d = edges,
  vertices = nodes %>% select(node_id, X, Y, NAMN1,ICC),
  directed = FALSE
)

g <- igraph::simplify( g,
  remove.multiple = TRUE,
  remove.loops    = TRUE)

  coord_dictionary <- nodes %>%
  # rename X→x and Y→y so the plot fn picks them up
  transmute(
    id   = as.character(node_id),
    x    = X,
    y    = Y,
    NAMN1,
    ICC
  )

stations_id <- coord_dictionary %>%
  filter(!is.na(NAMN1)) %>%
  pull(id)



  # NOw start simplyfying




non_protected <- V(g)[!(V(g)$name %in% get_protected_nodes(g,stations_id))]$name

condition <-any(degree(g,non_protected)==1)
while (condition){
    deg_one_nonprotected <- names(which(degree(g,non_protected)==1))
    g<- delete_vertices(g,deg_one_nonprotected)

    non_protected <- V(g)[!(V(g)$name %in% get_protected_nodes(g,stations_id))]$name
    condition <-any(degree(g,non_protected)==1)
}
graphs<-simplify_graph(g)
k= length(graphs)



  new_g <-collapse_contractible_Ys(g = graphs[[k]],stations_id = stations_id,maxdist = 3)

  final_gs <-simplify_graph(new_g)
t= length(final_gs)
final_g <- final_gs[[t]]


non_protected <- V(final_g)[!(V(final_g)$name %in% get_protected_nodes(final_g,stations_id))]$name

condition <-any(degree(final_g,non_protected)==1)
while (condition){
    deg_one_nonprotected <- names(which(degree(final_g,non_protected)==1))
    final_g<- delete_vertices(final_g,deg_one_nonprotected)

    non_protected <- V(final_g)[!(V(final_g)$name %in% get_protected_nodes(final_g,stations_id))]$name
    condition <-any(degree(final_g,non_protected)==1)
}


library(countrycode)   # iso-2/iso-3 ⇄ country names
## ── 1.  NODES TABLE ──────────────────────────────────────────────────────
nodes_df <- as_data_frame(final_g, what = "vertices") %>%
  transmute(
    nodeID        = name,             # keep igraph vertex name
    nodeLabel     = NAMN1,            # may be NA
    latitude      = Y,
    longitude     = X,
    iso2          = map_chr(ICC, ~ .x[[1]] %||% NA_character_),   # first code
    country_ISO3  = countrycode(iso2, "iso2c", "iso3c"),
    country_name  = countrycode(country_ISO3, "iso3c", "country.name")
  ) %>%
  select(nodeID, nodeLabel, latitude, longitude,
         country_name, country_ISO3)

## ── 2.  EDGES TABLE ──────────────────────────────────────────────────────
edges_df <- as_data_frame(final_g, what = "edges") %>%
  transmute(
    nodeID_from = from,     # already vertex names
    nodeID_to   = to
  )


## ── 3.  WRITE TO DISK ────────────────────────────────────────────────────
write.csv(nodes_df, "results/EE_nodes.csv",  row.names = FALSE, na = "")
write.csv(edges_df, "results/EE_edges.csv",  row.names = FALSE)


}