var nodes;
var arcs;
var node_scale;
var min_size;
var max_size;
var treestyle;
var reverse;
var highlight=-1;
var container;
var svgns = "http://www.w3.org/2000/svg";
var zoomTree;
var tree;
var cs_select;
var rescaled;
var from_value = [0,0,0,0];
var to_value = [1,1,1,1];
var min_value = [0,0,0,0];
var max_value = [1,1,1,1];
var pan = [-1,-1];
var selected_node;
var svgns = "http://www.w3.org/2000/svg";
var container;
var zoomTree;
var tree;
var cs_select;
var rescaled;
var highlight=-1;
var div_list=["properties"];

 //var slider = document.getElementById("myRange");
var textName;

 function downloadSVG() {
	document.getElementById("col1").children[0].appendChild(document.getElementById("grad0").children[0]);
	document.getElementById("col1").children[0].style.border='';
	document.getElementById("col1").children[0].children[1].setAttribute("x",30);
	document.getElementById("col1").children[0].children[1].setAttribute("y",50);
	var svg = container.cloneNode(true);
	svg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
	svg.setAttribute('onclick','');
	svg.setAttribute('ondblclick','');
	const blob = new Blob([svg.outerHTML.toString()]);
	const element = document.createElement("a");
	element.download = "JuDGEoutput.svg";
	element.href = window.URL.createObjectURL(blob);
	element.click();
	element.remove();
	document.getElementById("grad0").appendChild(document.getElementById("col1").children[0].children[1]);
	document.getElementById("col1").children[0].style.border='2px solid black';
}

	function rgb2hsl(r,g,b) {
		let v=Math.max(r,g,b), c=v-Math.min(r,g,b), f=(1-Math.abs(v+v-c-1));
		let h= c && ((v==r) ? (g-b)/c : ((v==g) ? 2+(b-r)/c : 4+(r-g)/c));
		return [60*(h<0?h+6:h), f ? c/f : 0, (v+v-c)/2];
	}

	function hsl2rgb(h,s,l) {
		let a= s*Math.min(l,1-l);
		let f= (n,k=(n+h/30)%12) => l - a*Math.max(Math.min(k-3,9-k,1),-1);
		return [f(0),f(8),f(4)];
	}

	function rgb2hex(rgb) {
		var r=rgb[0]*255;
		var g=rgb[1]*255;
		var b=rgb[2]*255;
		r=Math.round(r);
		g=Math.round(g);
		b=Math.round(b);
		return "#" + ((1 << 24) + (r << 16) + (g << 8) + b).toString(16).slice(1);
	}

	function convert2color(color_knots,value,scale_min,scale_max) {
		var temp=((value-scale_min)/(scale_max-scale_min)+0.005)/(1.01);
		if (reverse.checked) {
			temp=1-temp;
		}

		if (temp<0 || temp>1) {
			return color_knots.NA;
		}
		color_knots=color_knots.colors;
		for (var j=0;j<color_knots.length;j++) {
			if (temp<=color_knots[j]["weight"]) {
				var col;
				if (j==0) {
					col=color_knots[0]["color"];
				}
				else {
					var a=rgb2hsl(color_knots[j]["color"][0],color_knots[j]["color"][1],color_knots[j]["color"][2]);
					var b=rgb2hsl(color_knots[j-1]["color"][0],color_knots[j-1]["color"][1],color_knots[j-1]["color"][2]);

					var temp2=(temp-color_knots[j-1]["weight"])/(color_knots[j]["weight"]-color_knots[j-1]["weight"]);
					var h;
					if (a[0]-b[0]>180) {
						a[0]-=360;
						h=a[0]*temp2+b[0]*(1-temp2);
						if (h<0) {
							h+=360;
						}
					}
					else if (b[0]-a[0]>180) {
						b[0]-=360;
						h=a[0]*temp2+b[0]*(1-temp2);
						if (h<0) {
							h+=360;
						}
					}
					else {
						h=a[0]*temp2+b[0]*(1-temp2);
					}
					h=a[0]*temp2+b[0]*(1-temp2);
					col=hsl2rgb(h,a[1]*temp2+b[1]*(1-temp2),a[2]*temp2+b[2]*(1-temp2));
				}
				return col;
			}
		}

	}

	/*slider.oninput = function() {
		for (var i=0;i<8;i++) {
			data_keys[i]='';
		}
		highlight=-1;
		drawTemplate(this.value);
		show_details(selected)
		redrawFixed();
	 }*/

	var network=null;
	var selected=1;

  var svg1 = '<svg xmlns="http://www.w3.org/2000/svg" width="60" height="60" ondragover="onDragOver(event);" ondrop="onDrop(event);">content</svg>';

	var colors = ['#ED8181','#81DD81','#FDA04D','#9191FD'];
	var data_keys = ['','','','','','','',''];
	var x=30
	var y=30
	var r=28
	var el_index=-1
	var from_value = [0,0,0,0];
	var to_value = [1,1,1,1];
	var min_value = [0,0,0,0];
	var max_value = [1,1,1,1];
	var pan = [0,0];
	var expanded = false;

	function toggleCheckboxes() {
	  var checkboxes = document.getElementById("checkboxes");
	  if (!expanded) {
		checkboxes.style.display = "block";
		expanded = true;
	  } else {
		checkboxes.style.display = "none";
		expanded = false;
	  }
	}

	function confirmCheckboxes() {
	  var checkboxes = document.getElementById("checkboxes");
	  checkboxes.style.display = "none";
	  expanded = false;
		for (d in div_list) {
			c = document.getElementById("custom_chk_"+d);
			custom_div = document.getElementById(div_list[d]);
			if (c.checked) {
				custom_div.style.display = "block";
			}
			else {
				custom_div.style.display = "none";
			}
		}
	}

	function makeSegment(i,n,cx,cy,r_) {
		if ((n==1 && i==0) || i==-1) {
			segment = document.createElementNS(svgns, 'circle');
			segment.setAttributeNS(null, 'cx', cx);
			segment.setAttributeNS(null, 'cy', cy);
			segment.setAttributeNS(null, 'r', r_);
			return segment;
		}
		else if (i<n) {
			segment = document.createElementNS(svgns, 'path');
			var a=cx+r_*Math.cos(2*Math.PI/n*(i-0.5))
			var b=cy+r_*Math.sin(2*Math.PI/n*(i-0.5))

			var c=cx+r_*Math.cos(2*Math.PI/n*(i-1.5))
			var d=cy+r_*Math.sin(2*Math.PI/n*(i-1.5))

			var e = 'M' +	a.toString() + ',' +b.toString() +	'A' + r_.toString() +
				',' + r_.toString() +' 0 0,0 ' + c.toString() + ',' + d.toString() +
				' L ' + cx.toString() + ' ' + cy.toString() + ' z';

			segment.setAttributeNS(null, 'd', e);
			return segment;
		}
	}

	function makeSegment3(i,n,cx,cy,r_) {

		if ((n==1 && i==0) || (n==2 && i==0) || i==-1) {
			segment = document.createElementNS(svgns, 'circle');
			segment.setAttributeNS(null, 'cx', cx);
			segment.setAttributeNS(null, 'cy', cy);
			segment.setAttributeNS(null, 'r', r_);
			return segment;
		}
		else if (i<n-1) {
			segment = document.createElementNS(svgns, 'path');
			var a=cx+r_*Math.cos(2*Math.PI/(n-1)*(i-0.5))
			var b=cy+r_*Math.sin(2*Math.PI/(n-1)*(i-0.5))

			var c=cx+r_*Math.cos(2*Math.PI/(n-1)*(i-1.5))
			var d=cy+r_*Math.sin(2*Math.PI/(n-1)*(i-1.5))

			var e = 'M' +	a.toString() + ',' +b.toString() +	'A' + r_.toString() +
				',' + r_.toString() +' 0 0,0 ' + c.toString() + ',' + d.toString() +
				' L ' + cx.toString() + ' ' + cy.toString() + ' z';

			segment.setAttributeNS(null, 'd', e);
			return segment;
		}
		else {
			cr=r_/Math.sqrt(i+1.5);
			segment = document.createElementNS(svgns, 'circle');
			segment.setAttributeNS(null, 'cx', cx);
			segment.setAttributeNS(null, 'cy', cy);
			segment.setAttributeNS(null, 'r', cr);
			return segment;
		}
	}

	function makeSegment2(i, n, color) {
		if (color==null) {
			color='#EEEEEE';
		}

		if (n==1 && i==0) {
			var cr=r
			cr*=(n+1-i)/(n+1)
			el_index+=1;
			return '<circle id="a'+el_index.toString()+'" cx="'+x.toString()+'" cy="'+y.toString()+'" r="'+cr.toString()+'" '+
					'fill="'+color+'" stroke-width="2" stroke="#ffffff"/>';
		}
		else if (i<n) {
			var a=x+r*Math.cos(2*Math.PI/n*(i-0.5))
			var b=y+r*Math.sin(2*Math.PI/n*(i-0.5))

			var c=x+r*Math.cos(2*Math.PI/n*(i-1.5))
			var d=y+r*Math.sin(2*Math.PI/n*(i-1.5))
			el_index+=1;
			return '<path id="a'+el_index.toString()+'" d="M' +
				a.toString()+','+b.toString() +
				'A'+r.toString()+','+r.toString()+' 0 0,0 '+c.toString()+','+d.toString() +
				' L '+x.toString()+' '+y.toString()+' z" fill="' +
				color + '" stroke-width="2" stroke="#ffffff"/>';
		}
		else {
			var cr=r
			if (n>0) {
				cr/=Math.sqrt(i+1.5);
			}
			el_index+=1;
			return '<circle id="a'+el_index.toString()+'" cx="'+x.toString()+'" cy="'+y.toString()+'" r="'+cr.toString()+'" '+
					'fill="'+color+'" stroke-width="2" stroke="#ffffff"/>';
		}
	}

	function makeSegment4(i, n, color) {
		if (color==null) {
			color='#EEEEEE';
		}

		if (n==1 && i==0) {
			var cr=r
			cr*=(n+1-i)/(n+1)
			el_index+=1;
			return '<circle id="a'+el_index.toString()+'" cx="'+x.toString()+'" cy="'+y.toString()+'" r="'+cr.toString()+'" '+
					'fill="'+color+'" stroke-width="2" stroke="#ffffff"/>';
		}
		else {
			var a=x+r*Math.cos(2*Math.PI/n*(i-0.5))
			var b=y+r*Math.sin(2*Math.PI/n*(i-0.5))

			var c=x+r*Math.cos(2*Math.PI/n*(i-1.5))
			var d=y+r*Math.sin(2*Math.PI/n*(i-1.5))
			el_index+=1;
			return '<path id="a'+el_index.toString()+'" d="M' +
				a.toString()+','+b.toString() +
				'A'+r.toString()+','+r.toString()+' 0 0,0 '+c.toString()+','+d.toString() +
				' L '+x.toString()+' '+y.toString()+' z" fill="' +
				color + '" stroke-width="2" stroke="#ffffff"/>';
		}
	}

function preTitle(text) {
  const textcontainer = document.createElement("pre");
  textcontainer.innerText = text;
  return textcontainer;
}

function redrawFixed() {
	set_colors(true);
}

function redrawFixed2() {
	set_colors(true);
	pan[0]=container.children[0].transform.baseVal[0].matrix.e;
	pan[1]=container.children[0].transform.baseVal[0].matrix.f;
	show_details(highlight+1);
}

function updateSelect(j, select_data, select_key) {
	//slider_div.style.display="none";
	if (data_keys[j]==select_data[select_data.selectedIndex].text) {
		return;
	}

	data_keys[j]=select_data[select_data.selectedIndex].text;
	var option;
	select_key.innerHTML='';

	if (!(data_keys[j] in nodes[0].data)) {
		return;
	}

	if (Object.keys(nodes[0].data[data_keys[j]]).length>0) {
		for (key in nodes[0].data[data_keys[j]]) {
			option=document.createElement('option');
			option.text=key;
			select_key.add(option);
		}
		select_key.selectedIndex=0
	}
}

function gradscale(id,ct,_min,_max) {
	var svg2 = '<svg xmlns="http://www.w3.org/2000/svg" width="100%" height="40"><style> .small { font: 12px sans-serif; }</style>content</svg>';
	var stop='<stop offset="W%" style="stop-color:rgb(R%,G%,B%)"/>';
	var text='<text x="X" y="Y" class="small" text-anchor="ANCHOR">TEXT</text>'
	var line=' <line x1="X" y1="15" x2="X" y2="20" stroke="black" />'
	var temp='';

	for (var i=0;i<colorTemplates[ct].colors.length;i++) {
		temp+=stop.replace('W',colorTemplates[ct].colors[i].weight*100).replace('R',String(colorTemplates[ct].colors[i].color[0]*100)).replace('G',String(colorTemplates[ct].colors[i].color[1]*100)).replace('B',String(colorTemplates[ct].colors[i].color[2]*100));

		if (i!=colorTemplates[ct].colors.length-1) {
			var weight=(colorTemplates[ct].colors[i].weight+colorTemplates[ct].colors[i+1].weight)/2
			var rgb;
			if (reverse.checked) {
				rgb=convert2color(colorTemplates[ct],1-weight,0.0,1.0);
			}
			else {
				rgb=convert2color(colorTemplates[ct],weight,0.0,1.0);
			}
			temp+=stop.replace('W',weight*100).replace('R',String(rgb[0]*100)).replace('G',String(rgb[1]*100)).replace('B',String(rgb[2]*100));
		}
	}

	var joined

	var g_offset = document.getElementById("tree_content").childNodes[0].width.baseVal.value*0.6-16;
	var g_width = document.getElementById("tree_content").childNodes[0].width.baseVal.value*0.4;


	if (reverse.checked) {
		joined='<defs><linearGradient id="grad-ID" x1="100%" y1="0%" x2="0%" y2="0%">STOPS</linearGradient></defs><rect x="X" y="0" width="W" height="15" fill="url(#grad-ID)"/>'.replace('STOPS',temp).replace('ID',id).replace('ID',id).replace('X',g_offset.toString()).replace('W',g_width.toString());
	}
	else {
		joined='<defs><linearGradient id="grad-ID" x1="0%" y1="0%" x2="100%" y2="0%">STOPS</linearGradient></defs><rect x="X" y="0" width="W" height="15" fill="url(#grad-ID)"/>'.replace('STOPS',temp).replace('ID',id).replace('ID',id).replace('X',g_offset.toString()).replace('W',g_width.toString());
	}

	temp=''
	num_ticks=5
	for (var i=0;i<num_ticks;i++) {
		var a="middle";
		if (i==0) {
			//a="start";
			temp+=line.replace('X',(g_width/(num_ticks-1)*i+g_offset)+0.5).replace('X',(g_width/(num_ticks-1)*i+g_offset+0.5));
		}
		else if (i==num_ticks-1) {
			//a="end";
			temp+=line.replace('X',(g_width/(num_ticks-1)*i+g_offset)-0.5).replace('X',(g_width/(num_ticks-1)*i+g_offset-0.5));
		}
		else {
			temp+=line.replace('X',(g_width/(num_ticks-1)*i+g_offset)).replace('X',(g_width/(num_ticks-1)*i+g_offset));
		}

		temp+=text.replace('X',(g_width/(num_ticks-1)*i+g_offset)).replace('Y',30).replace('TEXT',formatNumber(_min+(_max-_min)*i/(num_ticks-1),5)).replace('ANCHOR',a);

	}

	document.getElementById('grad'+id.toString()).innerHTML = svg2.replace('content',joined+temp);
}

function addClickHandler(grad, i) {
	grad.addEventListener("click",function(){
		slider_div=document.getElementById("slider_div");
		if (slider_div.style.display=="none" || cs_select!=i) {
			cs_select=i;
			slider_div.style.display="block";
			$("#slider").ionRangeSlider({
				skin: "modern",
				type: "double",
				drag_interval: true,
				grid: true,
				onChange: function (data) {
					from_value[cs_select]=data.from;
					to_value[cs_select]=data.to;
					set_colors(false);
				},
			});
			s=Math.pow(10,Math.floor(Math.log10(Math.max(Math.abs(max_value[i]),Math.abs(min_value[i]))))-1)/2;

			$("#slider").data("ionRangeSlider").update({
				min: min_value[i],
				max: max_value[i],
				from: from_value[i],
				to: to_value[i],
				step: formatNumber(s,3),
			});
		}
		else {
			slider_div.style.display="none";
		}
	});
}

function drawTemplate(n) {
	var other = n%2;
	var n = Math.floor(n/2)+1;
	el_index=-1;
	/*
	var temp='<g transform="scale(1 1)">';

	if(other) {
		for (var i = 0; i<n; i++) {
			temp+=makeSegment2(i,n-1,colors[i]);
		}
	}
	else {
		for (var i = 0; i<n; i++) {
			temp+=makeSegment4(i,n,colors[i]);
		}
	}
	//document.getElementById("symbol").innerHTML = svg1.replace('content',temp).replace('"60"','"60"').replace('"60"','"60"');*/

	var element = document.getElementById("settings_content");
	element.innerHTML = '';

	treestyle = document.createElement('input');
    treestyle.type = 'checkbox';
    treestyle.id = 'chkLR';
    treestyle.checked = true;
	treestyle.addEventListener("change", fix);
	temp = document.createElement('label');
		temp.setAttributeNS(null, 'for','chkLR');
		temp.innerText='Tree style';
	element.appendChild(treestyle);
	element.appendChild(temp);

	reverse = document.createElement('input');
    reverse.type = 'checkbox';
    reverse.id = 'chkreverse';
    reverse.checked = false;
	reverse.addEventListener("change", fix);
	temp = document.createElement('label');
		temp.setAttributeNS(null, 'for','chkreverse');
		temp.innerText='Reverse colour scale';
	element.appendChild(reverse);
	element.appendChild(temp);

	textName = document.createElement('input');
		textName.type = 'checkbox';
		textName.id = 'chktextnode';
		textName.checked = true;
	textName.addEventListener("change", fix);
	temp = document.createElement('label');
		temp.setAttributeNS(null, 'for','chktextnode');
		temp.innerText='Show node names';
	element.appendChild(textName);
	element.appendChild(temp);

	for (var i = 0; i<n; i++) {
		var block = document.createElement("div");
		block.className="colorgroup";
		var selection = document.createElement("select");
		var selection2 = document.createElement("select");
		var selection3 = document.createElement("select");

		var col = document.createElement("span");
		var ref = document.createElement("span");
		var ref2 = document.createElement("span");
		var grad = document.createElement("span");

		col.style.backgroundColor=colors[i];
		col.className="colorgroup";
		col.innerText='   ';

		ref.className="colorgroup";
		ref2.className="colorgroup";
		selection.id="select"+i.toString();
		selection.class="select_join";
		selection2.id="select_data"+i.toString();
		selection3.id="select_key"+i.toString();
		grad.id="grad"+i.toString();
		grad.style.cursor="pointer"
		grad.style.padding="0px"
		addClickHandler(grad,i);
		var option;

		for (cs in colorTemplates) {
			option=document.createElement('option');
			option.text=cs;
			selection.add(option);
		}

		for (var key in scale) {
			option = document.createElement("option");
			option.text = key;
			selection2.add(option);
		}

		selection.selectedIndex = 0;
		selection.addEventListener("change", redrawFixed);
		selection2.addEventListener("change", redrawFixed2);
		selection2.selectedIndex = 0;
		selection3.addEventListener("change", redrawFixed2);
		selection3.selectedIndex = 0;
		updateSelect(i, selection2,selection3)


		block.appendChild(col);
		block.appendChild(selection);
		block.appendChild(ref);
		block.appendChild(selection2);
		block.appendChild(ref2);
		block.appendChild(selection3);
		//block.appendChild(grad);
		element.appendChild(block);
		element.appendChild(grad);
	}
}

function setup_positions(selected) {
	selected_node=selected;
	highlight=-1;
	var tree_div=document.getElementById("tree_content");
	var tree_frame=document.getElementById("tree_frame");

	container = document.createElementNS(svgns, "svg");
	container.setAttribute ("width", tree_frame.style.width);
	container.setAttribute ("height", String(Number(tree_frame.style.height.substr(0,tree_frame.style.height.length-2))-72)+"px");
	container.setAttributeNS(null, 'style', 'padding: 0px;' );
	container.addEventListener('wheel', evt => evt.preventDefault());
	container.setAttributeNS(null, 'onclick',"clickEvent(event,0,1)");
	container.setAttributeNS(null, 'ondblclick',"setup_positions(1)");

	if (tree_div.innerHTML!='') {
		tree_div.removeChild(tree_div.childNodes[0]);
		tree_div.insertBefore(container,tree_div.childNodes[0]);
	}
	else {
		tree_div.appendChild(container);
		tree_div.appendChild(document.getElementById("grad0"))
	}

	chooseTree();
	createNodes();

	subnodes=[];
	subedges=[];

	var minY=nodes[selected-1].Y;
	var maxY=nodes[selected-1].Y;
	var maxX=nodes[selected-1].X;

	if(selected==1) {
		subnodes=nodes;
		subedges=edges;
		for (n=0;n<nodes.length;n++) {
			if (nodes[n].X>maxX) {
				maxX=nodes[n].X;
			}
			if (nodes[n].Y>maxY) {
				maxY=nodes[n].Y;
			}
			else if (nodes[n].Y<minY) {
				minY=nodes[n].Y;
			}
		}
	}
	else {
		changed=true;
		ante=[selected];
		post=[selected];
		var breadth=0;
		var depth=nodes[nodes.length-1].level;
		while(changed) {
			changed=false;
			for (e=0;e<edges.length;e++) {
				if (ante.includes(edges[e].to)) {
					if (!ante.includes(edges[e].from)) {
						ante.push(edges[e].from);
						if (treestyle.checked) {
							nodes[edges[e].from-1].Y=nodes[selected-1].Y;
						}
						changed=true;
					}
				}
				if(post.includes(edges[e].from)) {
					if (!post.includes(edges[e].to)) {
						post.push(edges[e].to);
						if (treestyle.checked) {
							if (nodes[edges[e].to-1].Y>maxY) {
								maxY=nodes[edges[e].to-1].Y;
							}
							else if (nodes[edges[e].to-1].Y<minY) {
								minY=nodes[edges[e].to-1].Y;
							}
							if (nodes[edges[e].to-1].X>maxX) {
								maxX=nodes[edges[e].to-1].X;
							}
						}
						changed=true;
					}
				}
			}
		}

		for (n=0;n<nodes.length;n++) {
			if (ante.includes(nodes[n].id) || post.includes(nodes[n].id)) {
				subnodes.push(nodes[n]);
				if(nodes[n].level==depth) {
					breadth+=1;
				}
				for (e=0;e<edges.length;e++) {
					if (edges[e].from==nodes[n].id || edges[e].to==nodes[n].id) {
						if (!subedges.includes(edges[e])) {
							subedges.push(edges[e]);
						}
					}
				}
			}
		}
	}

	if (treestyle.checked) {
		var diff=maxY-minY;
		var sc=1.0;
		rescaled=false;
		if (diff < 30*max_size) {
			if (diff!=0) {
				sc=Math.pow(30*max_size/diff,0.75);
			}
			diff = 30*max_size;
			rescaled=true;
		}

		for (var i = 0; i < subnodes.length; i++) {
			subnodes[i].X=subnodes[i].X*diff/maxX;
			subnodes[i].Y=subnodes[i].Y*sc;
		}
	}

	for (var i = 0; i < subedges.length; i++) {
		if (subnodes.includes(nodes[subedges[i].from-1]) && subnodes.includes(nodes[subedges[i].to-1])) {
			line = document.createElementNS(svgns, 'line');
			line.setAttribute('x1',nodes[subedges[i].from-1].X.toString());
			line.setAttribute('y1',nodes[subedges[i].from-1].Y.toString());
			line.setAttribute('x2',nodes[subedges[i].to-1].X.toString());
			line.setAttribute('y2',nodes[subedges[i].to-1].Y.toString());
			line.setAttribute("stroke", "black");
			if (treestyle.checked) {
				line.setAttribute("stroke-width", 2.0);
			}
			else {
				line.setAttribute("stroke-width", 2*Math.pow(node_scale,nodes[subedges[i].from-1].level/2));
			}
			container.appendChild(line);
		}
	}

	for (var i = 0; i < subnodes.length; i++) {
		if (treestyle.checked) {
			subnodes[i].circle.setAttributeNS(null, 'x', subnodes[i].X-5*min_size);
			subnodes[i].circle.setAttributeNS(null, 'y', subnodes[i].Y-1.25*min_size);
		}
		else {
			subnodes[i].circle.setAttributeNS(null, 'x', subnodes[i].X-1.2*max_size*Math.pow(node_scale,nodes[i].level));
			subnodes[i].circle.setAttributeNS(null, 'y', subnodes[i].Y-1.2*max_size*Math.pow(node_scale,nodes[i].level));
		}
		container.appendChild(subnodes[i].circle);
	}

	redrawFixed();

	delete zoomTree
	zoomTree = svgPanZoom(container, {zoomScaleSensitivity: 0.3,preventMouseEventsDefault: false, dblClickZoomEnabled: false});
	zoomTree.zoomBy(0.88)
	pan[0]=-1;
	pan[1]=-1;
	show_details(0);
}

function clickEvent(event,i,j) {
	if(event.ctrlKey) {
		setup_positions(j);
	}
	else {
		show_details(i);
	}
}

function createNodes()
{
	for (var i = 0; i < nodes.length; i++) {
		delete nodes[i].circle;
		nodes[i].circle = document.createElementNS(svgns, 'svg');
		if (treestyle.checked) {
			nodes[i].circle.setAttributeNS(null, 'width', 10*min_size);
			nodes[i].circle.setAttributeNS(null, 'height', 3.5*min_size);
		}
		else {
			nodes[i].circle.setAttributeNS(null, 'width', 2.4*max_size*Math.pow(node_scale,nodes[i].level));
			nodes[i].circle.setAttributeNS(null, 'height', 2.4*max_size*Math.pow(node_scale,nodes[i].level));
		}
		nodes[i].circle.setAttributeNS(null, 'cursor','pointer');
		nodes[i].circle.setAttributeNS(null, 'ondblclick',"event.stopPropagation();setup_positions("+(i+1).toString()+")");
		nodes[i].circle.setAttributeNS(null, 'onclick',"event.stopPropagation();clickEvent(event,"+(i+1).toString()+","+(i+1).toString()+")");
	}
}
  function set_colors(reset) {
	//var other=slider.value % 2;
	var other=0;
	var n=1;
	//var n = Math.floor(slider.value/2)+1;
	var element;
	var ct;
	var min_ = {};
	var max_ = {};

	if (reset) {
		for (var j=0;j<n;j++) {
			select_data=document.getElementById("select_data"+j.toString());
			select_key=document.getElementById("select_key"+j.toString());
			updateSelect(j, select_data, select_key)
			if (data_keys[j]!='') {
				element = document.getElementById("select"+j.toString());
				ct=element.options[element.selectedIndex].text;
				if (ct in max_) {
					if (scale[data_keys[j]].max>max_[ct]) {
						max_[ct]=scale[data_keys[j]].max;
					}
					if (scale[data_keys[j]].min<min_[ct]) {
						min_[ct]=scale[data_keys[j]].min;
					}
				}
				else {
					max_[ct]=scale[data_keys[j]].max;
					min_[ct]=scale[data_keys[j]].min;
				}
			}

		}
		for (var j=0;j<n;j++) {
			if (data_keys[j]!='') {
				element = document.getElementById("select"+j.toString());
				ct=element.options[element.selectedIndex].text;
				max_value[j]=parseFloat(formatNumber(max_[ct],2));
				if (max_value[j]<max_[ct]) {
					max_value[j]=parseFloat(formatNumber(max_[ct]+Math.pow(10,Math.floor(Math.log10(Math.abs(max_[ct])))-1)/2,2));
				}
				min_value[j]=parseFloat(formatNumber(min_[ct],2));
				if (min_value[j]>min_[ct]) {
					min_value[j]=parseFloat(formatNumber(min_[ct]-Math.pow(10,Math.floor(Math.log10(Math.abs(min_[ct])))-1)/2,2));
				}

				from_value[j]=min_value[j];
				to_value[j]=max_value[j];
			}
		}
	}

	for (var j=0;j<n;j++) {
		element = document.getElementById("select"+j.toString());
		ct=element.options[element.selectedIndex].text;
		gradscale(j,ct,from_value[j],to_value[j]);
	}

	for (var i = 0; i < nodes.length; i++) {
		var temp='';
		var dat_col=['#EEEEEE','#EEEEEE','#EEEEEE','#EEEEEE','#EEEEEE','#EEEEEE','#EEEEEE','#EEEEEE'];
		nodes[i].title='Node '+nodes[i].label;

		while (nodes[i].circle.lastChild) {
		nodes[i].circle.removeChild(nodes[i].circle.lastChild);
		}

		for (var j=-1;j<n;j++) {
			var subtext='';
			element = document.getElementById("select"+j.toString());
			if(j>=0) {
				ct=element.options[element.selectedIndex].text;
				dat_col[j]=rgb2hex(colorTemplates[ct].NA);
				if (data_keys[j]!='') {
					if (data_keys[j] in nodes[i].data) {
						nodes[i].title+='\n';
						if (Object.keys(nodes[i].data[data_keys[j]]).length==0) {
							dat_col[j]=rgb2hex(convert2color(colorTemplates[ct],nodes[i].data[data_keys[j]],from_value[j],to_value[j]));
							if (nodes[i].data[data_keys[j]]!=null) {
								subtext=formatNumber(nodes[i].data[data_keys[j]],6);
								nodes[i].title+=data_keys[j]+': '+subtext;
							}
						}
						else {
							select_key=document.getElementById("select_key"+j.toString());
							key=select_key[select_key.selectedIndex].text;
							dat_col[j]=rgb2hex(convert2color(colorTemplates[ct],nodes[i].data[data_keys[j]][key],from_value[j],to_value[j]));
							if (nodes[i].data[data_keys[j]][key]!=null) {
								subtext=formatNumber(nodes[i].data[data_keys[j]][key],6);
								nodes[i].title+=data_keys[j]+'('+key+'): '+subtext;
							}
						}
					}
				}
			}
			var circle;// = document.createElementNS(svgns, 'circle');
			sw=0
			var cx;
			var cy;
			var r_;
			if (treestyle.checked) {
				cx=5*min_size;
				cy=1.25*min_size;
				//circle.setAttributeNS(null, 'cx', 5*min_size);
				//circle.setAttributeNS(null, 'cy', 1.25*min_size);
				r_ = min_size;
				//circle.setAttributeNS(null, 'r', r_);
				sw = 2
			}
			else {
				cx=1.2*max_size*Math.pow(node_scale,nodes[i].level);
				cy=1.2*max_size*Math.pow(node_scale,nodes[i].level)
				//circle.setAttributeNS(null, 'cx', 1.1*max_size*Math.pow(node_scale,nodes[i].level));
				//circle.setAttributeNS(null, 'cy', 1.1*max_size*Math.pow(node_scale,nodes[i].level));
				r_ = max_size*Math.pow(node_scale,nodes[i].level);
				//circle.setAttributeNS(null, 'r', r_);
				sw=3*Math.pow(node_scale,nodes[i].level/2)
			}

			if(other)
				circle=makeSegment3(j,n,cx,cy,r_);
			else {
				circle=makeSegment(j,n,cx,cy,r_);
			}

			if (j==-1) {
				circle.setAttributeNS(null, 'style', 'fill: black; stroke: black; stroke-width: '+(sw/2).toString()+'px;' );
			}
			else {
				circle.setAttributeNS(null, 'style', 'fill: '+dat_col[j]+'; stroke: white; stroke-width: '+sw.toString()+'px;' );
			}
			nodes[i].circle.appendChild(circle);
			if (j==0 && treestyle.checked) {
				text=document.createElementNS(svgns, 'text');
				if (nodes[i].leaf && !rescaled) {
					text.setAttributeNS(null, 'x', 6.2*min_size);
					text.setAttributeNS(null, 'y', 1.35*min_size);
					text.setAttributeNS(null, 'text-anchor', 'left');
				}
				else {
					text.setAttributeNS(null, 'x', 5*min_size);
					text.setAttributeNS(null, 'y', 3.1*min_size);
					text.setAttributeNS(null, 'text-anchor', 'middle');
				}
				text.setAttributeNS(null, 'font-size', 24);
				text.setAttributeNS(null, 'dominant-baseline', 'middle');
				if(textName.checked) {
					text.appendChild(document.createTextNode(nodes[i].label));
				}
				else {
					text.appendChild(document.createTextNode(subtext));
				}

				nodes[i].circle.appendChild(text);
			}
		}


		var title = document.createElementNS(svgns,"title");
		title.textContent = nodes[i].title;
		nodes[i].circle.appendChild(title);

	}
}

function formatNumber(value,n) {
	value+=0.1
	value-=0.1
	sf=getSignificantDigitCount(value);
	if (sf>=n) {
		return value.toPrecision(n);
	}
	else {
		return value.toPrecision(sf);
	}
}

function getSignificantDigitCount(n) {
	n2 = Math.abs(n).toFixed(15).replace(".", ""); //remove decimal and make positive
	if (n2 == 0) return 1;
	while (n2 != 0 && n2 % 10 == 0) n2 /= 10; //kill the 0s at the end of n

	return Math.max(Math.floor(Math.log10(Math.abs(n))),Math.floor(Math.log10(n2))) + 1; //get number of digits
}

function show_details(selected) {
	if(selected>0 || (pan[0]==-1 && pan[1]==-1)) {
		pan[0]=container.children[0].transform.baseVal[0].matrix.e;
		pan[1]=container.children[0].transform.baseVal[0].matrix.f;
	}

	if(container.children[0].transform.baseVal[0].matrix.e!=pan[0] || container.children[0].transform.baseVal[0].matrix.f!=pan[1]) {
		pan[0]=container.children[0].transform.baseVal[0].matrix.e;
		pan[1]=container.children[0].transform.baseVal[0].matrix.f;
		return;
	}

	if(highlight>=0) {
		nodes[highlight].circle.children[0].style.stroke="white";
		temp=nodes[highlight].circle.children[0].style.strokeWidth;
		width=parseFloat(temp.substr(0,temp.length-2));
		nodes[highlight].circle.children[0].style.strokeWidth="0px";
	}
	send_message(selected.toString()+";"+document.getElementById("select_data0").value+";"+document.getElementById("select_key0").value,"default");

	if(selected>=1) {
		highlight=selected-1;
		nodes[highlight].circle.children[0].style.stroke="black"
		temp=nodes[highlight].circle.children[0].style.strokeWidth;
		width=parseFloat(temp.substr(0,temp.length-2));
		nodes[highlight].circle.children[0].style.strokeWidth="7px";
	}
	else {
		highlight=-1
	}

	if (typeof node_select === 'function') {
		node_select(selected);
	}
	/*confirmCheckboxes();

	}
	else {
		highlight=-1
		for (d in div_list) {
			custom_div = document.getElementById(div_list[d]);
			custom_div.style.display = "none";
		}
	}*/
}

function chooseTree() {
	if (treestyle.checked) {
		for (var i = 0; i < nodes.length; i++) {
			nodes[i].X=nodes[i].posX2;
			nodes[i].Y=nodes[i].posY2;
		}
	}
	else {
		for (var i = 0; i < nodes.length; i++) {
			nodes[i].X=nodes[i].posX;
			nodes[i].Y=nodes[i].posY;
		}
	}
}

function fix() {
	setup_positions(selected_node);
}

function create_custom_dropdown() {
	cb = document.getElementById("list_custom")
	for(d in div_list) {
		var l = document.createElement("label");
		var c = document.createElement("input");
		var s = document.createElement("span");
		c.setAttribute("type", "checkbox");
		c.id = "custom_chk_"+d;
		c.checked=true;
		l.appendChild(c);
		s.setAttribute("class", "chklabel");
		s.innerHTML = div_list[d];
		l.appendChild(s);
		cb.appendChild(l);
	}
}

function load_jsonview(selected,div) {

	if (selected==0) {
		div.innerHTML="No node selected";
		div.style.padding="10px";
		div.parentElement.style.overflow="auto";
		return;
	}

	var temp=JSON.parse(JSON.stringify(nodes[selected-1]));
	delete temp["id"];
	delete temp["level"];
	delete temp["posX"];
	delete temp["posY"];
	delete temp["posX2"];
	delete temp["posY2"];
	delete temp["circle"];
	delete temp["title"];
	delete temp["X"];
	delete temp["Y"];
	delete temp["leaf"];

	temp["Name"]=temp["label"];
	delete temp["label"];

	for ( obj in temp["data"]){
		temp[obj]=temp["data"][obj];
	}
	delete temp["data"];

	tree = JsonView.createTree(temp);

	// render tree into dom element
	div.innerHTML="";
	JsonView.render(tree, div);
	JsonView.expandChildren(tree);
	//JsonView.toggleNode(tree);
}
