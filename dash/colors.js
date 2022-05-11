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
	var temp=(value-scale_min)/(scale_max-scale_min)
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

var colorTemplates = {blackredgreenblue:{NA:[0.9,0.9,0.9],colors:[{weight:0.0,color:[0.0,0.0,0.0]},{weight:0.333,color:[0.85,0.0,0.0]},{weight:0.667,color:[0.0,0.8,0.0]},{weight:1.0,color:[0.0,0.0,0.85]}]},
blackgreenblue:{NA:[0.9,0.9,0.9],colors:[{weight:0.0,color:[0.0,0.0,0.0]},{weight:0.5,color:[0.0,0.8,0.0]},{weight:1.0,color:[0.0,0.0,0.85]}]},
black:{NA:[0.9,0.9,0.9],colors:[{weight:0.0,color:[0.0,0.0,0.0]},{weight:1.0,color:[0.0,0.0,0.0]}]},
blackwhite:{NA:[0.7,0.7,0.9],colors:[{weight:0.0,color:[0.0,0.0,0.0]},{weight:1.0,color:[1.0,1.0,1.0]}]},
greenblue:{NA:[0.9,0.9,0.9],colors:[{weight:0.0,color:[0.0,0.7,0.5]},{weight:1.0,color:[0.0,0.0,0.9]}]},
blackred:{NA:[0.9,0.9,0.9],colors:[{weight:0.0,color:[0.0,0.0,0.0]},{weight:1.0,color:[0.95,0.0,0.0]}]},
bgyr:{NA:[0.9,0.9,0.9],colors:[{weight:0.0,color:[0.0,0.0,0.95]},{weight:0.5,color:[0.0,0.95,0.0]},{weight:1.0,color:[0.95,0.0,0.0]}]},
bpbcgyr:{NA:[0.9,0.9,0.9],colors:[{weight:0.0,color:[0.0,0.0,0.0]},{weight:0.1667,color:[0.95,0.0,0.95]},{weight:0.3333,color:[0.0,0.0,0.95]},{weight:0.5,color:[0.0,0.95,0.95]},{weight:0.6667,color:[0.0,0.95,0.0]},{weight:0.8333,color:[0.95,0.95,0.0]},{weight:1.0,color:[0.95,0.0,0.0]}]}
};
