function js_cb1(data) {
	if (data) {
		send_message("The selected node is a leaf node","custom");
	}
	else {
		send_message("The selected node is not a leaf node","custom");
	}
}

function js_cb2(data) {
	document.getElementById("block1_content").style.color=data;
}
