function node_select(selected) {
	send_message([selected,document.getElementById("select_data0").value,document.getElementById("select_key0").value],"default");
}

function js_cb1(data) {
	if (data==null) {
		send_message("There is no selection","custom");
	}
	else if (data==true) {
		send_message("The selected node is a leaf node","custom");
	}
	else {
		send_message("The selected node is not a leaf node","custom");
	}
}

function js_cb2(data) {
	document.getElementById("block3").style.color=data;
}
