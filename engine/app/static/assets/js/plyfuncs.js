var d3 = Plotly.d3;

//
// 2D Scatter
//

var myPlot = document.getElementById('myDivPlot');

var xmin = +1E+8;
var ymin = +1E+8;
var xmax = -1E+8;
var ymax = -1E+8;

for ( trace_id in DATA ) {
	var x = DATA[trace_id].x;
	var y = DATA[trace_id].y;
    if (Math.min.apply(null,x) < xmin){
    	xmin = Math.min.apply(null,x);
    }
    if (Math.max.apply(null,x) > xmax) {
    	xmax = Math.max.apply(null,x);
    }
	if (Math.min.apply(null,y) < ymin){
    	ymin = Math.min.apply(null,y);
    }
    if (Math.max.apply(null,y) > ymax) {
    	ymax = Math.max.apply(null,y);
    }
}

for ( trace_id in DATA ) {

	var trace_scatter = {
		name: DATA[trace_id].name,
		text: DATA[trace_id].label,
		x: DATA[trace_id].x,
		y: DATA[trace_id].y,
		ids: DATA[trace_id].id,
		type: 'scatter',
		mode: 'markers',
		marker: {size:10, opacity: 1.0}
	};

	var trace_density = {
		x: DATA[trace_id].x,
		y: DATA[trace_id].y,
		name: 'density',
		ncontours: 20,
		colorscale: 'Hot',
		reversescale: true,
		showscale: false,
		type: 'histogram2dcontour'
	};

	var layout = {
		title: 'QCB Analysis',
		xaxis: {title: LABX, range:[xmin,xmax], titlefont: {size: 24}, tickfont: {size: 18}},
		yaxis: {title: LABY, range:[ymin,ymax], titlefont: {size: 24}, tickfont: {size: 18}},
		hovermode: 'closest',
		dragmode: 'lasso'
		// shapes: [{
		// 	type: 'line',
		// 	x0: 0,
		// 	y0: 0,
		// 	x1: 2.0E6,
		// 	y1: 0.4E6,
		// 	line: {color:'rgb(128, 0, 128)', width:4, dash:'dot'}
		// }]
	};

	Plotly.plot(myPlot, [trace_scatter], layout);

}

//
// Output Table
//



var myTable = document.getElementById('myDivTable');

var table_x = [];
var table_y = [];
var table_i = [];

for ( trace_id in DATA ) {
	for (pt = 0; pt < DATA[trace_id].name.length; pt++) {
		table_x.push(DATA[trace_id].x[pt]);
		table_y.push(DATA[trace_id].y[pt]);
		table_i.push(DATA[trace_id].id[pt]);
	}
}

var trace = {
	type: 'table',
	header: {values: [["id"],["x"],["y"]]},
	cells: {values: [table_i,table_x, table_y]},
	columnwidth: [400,200,200],
};

Plotly.plot(myTable, [trace], {margin: {l: 5, r: 5, b: 5, t: 5, pad: 0}});

myPlot.on('plotly_selected', function(data){
		if (data != null) { // Does not go in when clicked
				var  cid = [];
				var xsel = [];
				var ysel = [];
				data.points.forEach(function(pt) {
						cid.push(pt.id);
						xsel.push(pt.x);
						ysel.push(pt.y);
				});
				Plotly.restyle(myTable, {
						cells: {values: [cid,xsel,ysel]}
				}, [0]);
		}
});

function updateImgs() {
		for (j in myTable.data[0].cells.values[0]) {
				var img_name = myTable.data[0].cells.values[0][j]
				var img = document.createElement("IMG");
				// img.setAttribute("width", "128");
				img.setAttribute("src","static/imgs/"+img_name+".jpg?" + new Date().getTime());
				//img.setAttribute("src","https://cran.r-project.org/web/packages/plotly/readme/man/figures/plotly.png");
				img.setAttribute("title",img_name);
				document.getElementById('myDivImgs').appendChild(img);
		}
}

function clearImgs() {
		var list = document.getElementById("myDivImgs");
		while (list.firstChild) {
				list.removeChild(list.firstChild);
		}
		Plotly.restyle(myDivPlot, {selectedpoints: [null]});
}
