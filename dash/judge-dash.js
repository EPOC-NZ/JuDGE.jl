function initialize(data) {
	nodes = data.nodes;
	edges = data.arcs;
	node_scale = data.node_scale;
	min_size = data.min_size;
	max_size = data.max_size;
	scale = data.scale;

	drawTemplate(0);
	setup_positions(1);
}
