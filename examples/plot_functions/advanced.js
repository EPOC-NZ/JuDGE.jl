function node_select(selected) {
	load_jsonview(selected,document.getElementById("block6"));

	send_message([selected,document.getElementById("select_data0").value,document.getElementById("select_key0").value],"default");
}
