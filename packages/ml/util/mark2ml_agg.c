main() {

  int j, k, agg_0, node, list_node[300], Nlist = 0, current_agg = -1;
  while ( scanf("%d%d",&agg_0, &node) == 2 ) {
    if (agg_0 != current_agg) {
      for (j = 0; j < Nlist; j++) printf("%d %d\n",list_node[j],current_agg);
      current_agg = agg_0;
      Nlist = 1;
      list_node[0] = node;
    }
    else {
      for (k = 0; k < Nlist; k++) {
	if (list_node[k] > node) break;
      }
      for (j = Nlist-1; j >= k; j--) 
	list_node[j+1] = list_node[j];
      list_node[k] = node;
      Nlist++;
    }
  }
  for (j = 0; j < Nlist; j++) printf("%d %d\n",list_node[j],current_agg);
}
