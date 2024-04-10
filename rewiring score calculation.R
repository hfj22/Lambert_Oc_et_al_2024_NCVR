library(igraph)

network<-read.table("nonrsp_prp_GRNs.csv", sep=",")
colnames(network)<-network[1,]
network<-network[2:nrow(network),]

network_prp<-network[,-c(3)]
network_nonrsp<-network[,1:3]

n_prp<-graph_from_data_frame(network_prp, directed = TRUE)
n_nonrsp<-graph_from_data_frame(network_nonrsp, directed = TRUE)

adj_prp<-as_adjacency_matrix(n_prp,attr= "prp_coef_mean",sparse=FALSE)
adj_nonrsp<-as_adjacency_matrix(n_nonrsp,attr= "nonrsp_coef_mean",sparse=FALSE)

adj_prp_num<-matrix(as.numeric(adj_prp), ncol=ncol(adj_prp), dimnames = list(rownames(adj_prp), colnames(adj_prp)))
adj_nonrsp_num<-matrix(as.numeric(adj_nonrsp), ncol=ncol(adj_nonrsp), dimnames = list(rownames(adj_nonrsp), colnames(adj_nonrsp)))

for (i in 1:nrow(adj_prp_num)) {
  for (j in 1:ncol(adj_prp_num)) {
    if (is.na(adj_prp_num[i,j])==TRUE) {
      adj_prp_num[i,j] = 0
    }
  }
}

for (i in 1:nrow(adj_nonrsp_num)) {
  for (j in 1:ncol(adj_nonrsp_num)) {
    if (is.na(adj_nonrsp_num[i,j])==TRUE) {
      adj_nonrsp_num[i,j] = 0
    }
  }
}

delta_adj<-adj_prp_num-adj_nonrsp_num
nodes_score<-as.data.frame(rowSums(abs(delta_adj)))
colnames(nodes_score)<-"out_score"
nodes_score$in_score<-colSums(abs(delta_adj))
nodes_score$total_score<-nodes_score[,1]+nodes_score[,2]

write.table(nodes_score, "rewiring scores.csv")
