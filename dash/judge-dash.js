 var cb_sleep={};

function callback_monitors(cb_names) {
	for (var i=0;i<cb_names.length;i++) {
		callback_monitor(cb_names[i]);
	}
}

function callback_monitor(cb_name) {
	cb_sleep[cb_name]=false;
	var target = document.getElementById('cb-'+cb_name+'_content');
	var observer = new MutationObserver(function(mutations) {
		var data = JSON.parse(sessionStorage["data-"+cb_name]);

		cb_sleep[cb_name]=!cb_sleep[cb_name];

		if (cb_sleep[cb_name]) {
			return;
		}
		callback_controller(cb_name,data);
		sessionStorage["data-"+cb_name]="";
	});

	observer.observe(target, {
		attributes:    true,
		childList:     false,
		characterData: false
	});
}

$( document ).ready(function() {
   setTimeout(function() {
	var tree_data = JSON.parse(sessionStorage["tree_data"]);
	sessionStorage["tree_data"]="";
	sessionStorage["data_default"]="";
	nodes = tree_data.nodes;
	edges = tree_data.arcs;
	node_scale = tree_data.node_scale;
	min_size = tree_data.min_size;
	max_size = tree_data.max_size;
	scale = tree_data.scale;
	cbs = tree_data.callbacks;
	drawTemplate(0);

	setup_positions(1);

	if (cbs.length!=0) {
		callback_monitors(cbs);
	}

}, 500);
});

function send_message(message,target) {
	temp=store.getState().paths.strs["msg-"+target];
	i=temp[temp.length-1];

	store.getState().layout.props.children.props.children[i].props.data=message.toString();
	sessionStorage["msg-"+target]=message.toString();
	document.getElementById("cb-"+target+"-button").click();
}
