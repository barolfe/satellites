<?php
	$TLEData = file_get_contents('http://celestrak.com/NORAD/elements/stations.txt');
	$TLEData .= file_get_contents('http://www.celestrak.com/NORAD/elements/visual.txt');
	$TLEData .= file_get_contents('http://www.celestrak.com/NORAD/elements/molniya.txt');
	$TLEData .= file_get_contents('http://www.celestrak.com/NORAD/elements/sbas.txt');
	$TLEData .= file_get_contents('https://celestrak.com/NORAD/elements/geodetic.txt');
	$TLEData .= file_get_contents('http://www.celestrak.com/NORAD/elements/geo.txt');
	$TLEData .= file_get_contents('http://www.celestrak.com/NORAD/elements/weather.txt');
	$TLEData .= file_get_contents('http://www.celestrak.com/NORAD/elements/x-comm.txt');
	$TLEData .= file_get_contents('http://www.celestrak.com/NORAD/elements/tle-new.txt');
	$TLEData .= file_get_contents('http://www.celestrak.com/NORAD/elements/resource.txt');

?>

<!Doctype HTML>
<html lang="en">

<head>
<meta charset="UTF-8">
<!-- Some phones will try to autoformat telephone numbers, which interferes with the two-line element data.
	This is a way to turn that off. -->
<meta name="format-detection" content="telephone=no"> 

<title>BryanRolfe.com - Look up!</title>

<style type="text/css">
	body {
		background:#000000; 
		font-family: sans-serif; 
		color: rgb(200, 200, 200); 
		font-size: 14px; 
		transition: background 1s;}
	a {color: inherit;}
	p {font-family: sans-serif; font-size: 14px; margin: 4px;}
	table {text-align: left; width: 100%;}

    select option {
    	color: #000000;
    }

	.top { position: relative; width: 800px; height:auto; left: 50%; margin-left: -400px}
	.blog {position: relative; width: 800px; left:50%; margin-left:-400px; padding-top:20px; padding-bottom: 20px; text-align: left}
	.blog p { padding-top: 10px;}
	#objectInfo {position: relative; width: 700px; left:50%; margin-left:-350px; display:block;}
	.minheight{
		min-height: 14px;
		clear: both;
	}
	.left {float:left;}
	.right {float:right; }
	.coords {background: transparent; color: inherit; font-size:inherit; border: 0px; width: 80px;}

	#activeSatellites {
	}
	#activeSatellites a {
		margin-left: 5px;
	}

	#earthContainer {
		width: 800px;
		position: relative;
	}
	#earthInset {
		float: right;
		bottom: 0px;
		left: 15px;
		position: absolute;
	}
	#close {
		padding: 5px;
		text-decoration: none;
		right: 0px;
	}
	#canvasHolder {
		border: 1px dashed;
		border-color: rgba(255, 255, 255, 0.3); 
		/* background: rgba(0, 0, 0, 0.5); */
	}

	#pointer {
		position: absolute;
		margin-left:-130px;
		padding: 3px;
	}
	select {
		color: inherit;
		font-size: inherit;
		background: transparent;
		padding: 3px;
		margin: 0;
		display: inline-block;
	    border: none;
	    outline: none;
	    -webkit-appearance:none;
	    -moz-appearance:none;
	    appearance:none;
		cursor: pointer;
	}

</style>


</head>

<body style="text-align:center">
<div id="TLEData" style="display: none;"><?php echo $TLEData; ?></div>

<div id="earthContainer" style="margin-left:auto; margin-right:auto">
	<a href="#" style="text-decoration: none;" onclick="toggle_color(); init(); satellite(); mercator();">Change Color</a></p>
	<!-- <canvas id="satelliteCanvas" width="400" height="400"></canvas> -->
	<svg id="earth">

	</svg>
	<div id="earthInset">
		<a id="close" href="#" style="float:left">X</a>
		<div id="canvasHolder">
		</div>
	</div>
</div>



<div class = "top">
	The code for this web app is open source and can be found at: <a href="https://github.com/barolfe/satellites">https://github.com/barolfe/satellites</a><br>
	Where are you?
	<div class="minheight">
		Longitude: <input id = "uLongi" class="coords" type="text" value="0.00"></input>
		Latitude: <input id = "uLat" class="coords" type="text" value="0.00"></input>
	</div>
	Time Future (debug mode): <input id = "timeForward" class="coords" type="text" value="0.00"></input>
	<div id="activeSatellites" class="minheight">
		<!-- This area is populated by user-selected active satellites -->
	</div>
	<div class="minheight">
		<div id="pointer">Pick a satellite --></div>
		<select id="selectSatellites" style="float:left">
			<option value="-1">Pick Satellite</option>
		</select>
		 <p class="left">[ <a href="#" onclick="toggle_visibility('objectInfo')">show/hide info</a> ]</p>
	</div>
	<div id="objectInfo" style="display:block">
		<div class="minheight">
			<table style="width:100%; font-size: 14px">
				<tr>
					<td>Inclination</td> <td> Argument of Perigee </td> <td>Ascending Node</td> <td>Mean Anomoly</td> <td>Mean Motion</td> <td>Eccentricity</td>
				</tr>
				<tr>
					<td id="inclination"></td> <td id="aoPerigee"></td> <td id="raAscending"></td> <td id="meanAnomoly"></td> <td id="meanMotion"></td> <td id="eccentricity"></td>
				</tr>
			</table>
		</div>
		<div class="minheight"> </div>
		<div class="minheight">
			<p class="left"> Last NORAD Sync:</p> <p id="epochTime" class="right"></p>
		</div>
		<div class="minheight">
			<p class="left"> Current Time:</p> <p id="currentTime" class="right"></p>
		</div>
		<div class="minheight">
		</div>
		<div class="minheight">
			<p class="left"> Eccentric Anomoly:</p> <p id="eccentricAnomoly" class="right"></p>
		</div>
		<div class="minheight">
			<p class="left"> Velocity:</p> <p id="velocity" class="right"></p>
		</div>
		<div class="minheight">
			<p class="left"> Altitude:</p> <p id="altitude" class="right"></p>
		</div>
	</div>

</div>


<div class="blog">
<p style="font-weight: bold;">A plea to look up</p>
<p>
Of all the good that came with the industrial revolution and the invention of the electric light bulb, we lost something that may not be apparent -- the night sky. In most developed communities, a glance towards the heavens reveals only a smattering of stars -- the ones bright enough to pierce through the glare that gives us our night time productivity, safety, and comfort. 
</p>
<p>
It is perhaps ironic that the same forces which hid the stars also spawned the race to space, ultimately lighting the night sky with our own stars. But these new stars are not luminous balls of plasma, rather they are satellites. One such satellite deserves special attention, not only is it visible from the most light-polluted cities, it is also the only place humanity continually calls home beyond Earth: the International Space Station.
</p>
<p>
Given a dark location far from city lights, you can easily see the Milky Way, including its vast nebula and a large region of dust and gas known as the Great Rift. In fact, you can even see our galactic neighbor, Andromeda, a spiral galaxy a mere 2.5 million light years away which is currently on a collision course with the Milky Way, our own galaxy. If you throw binoculars into the mix, other galaxies like M51 and M83 (not to be confused with the French electronic band) also become visible. It is a humbling experience and undeniably the greatest show on, and off, Earth. 
</p>
<p>
Unfortunately, finding dark skies these days is a difficult task. A quick peak at a light pollution map reveals that such skies are becoming increasingly elusive. There are things that can be done to mitigate light pollution, and there are movements dedicated to that cause, but I am not writing to address those issues. I am writing to turn your head skywards and spread awareness of something that everyone can see.
</p>
<p>
Traveling at over <span id="ISSspeed"></span> at an altitude of <span id="ISSalt"></span> above the Earth, six humans skim the atmosphere aboard a one million pound structure larger than a football field. Aboard the International Space Station, these fortunate few experience a day that is a mere 92 minutes long with 15 sunrises and sunsets per Earth-day. They arguably have the best view of the night sky and the only view of Earth from space. 
</p>
<p>
For astronauts like Chris Hadfield, the ISS becomes a home, a laboratory, and the perfect backdrop to cover David Bowie’s <a href="http://www.youtube.com/watch?v=KaOC9danxNo">Space Oddity</a>. But for the rest of us, it’s an alien environment inaccessible due to the extreme difficulty, cost, and dangers of orbital launches. The good news is that everyone can be a part of the ISS experience.
</p>
<p>
Although you likely won’t have the chance to orbit our planet, you can still see the ISS without the need for a telescope or even binoculars. If the timing is right, the ISS rivals even the brightest stars in the night sky. If you can see Vega, one of the brightest stars in the sky, or Jupiter, you can see the ISS. That’s not terribly surprising when you consider that it’s a football field-sized reflective structure bouncing sunlight to our eyes. 
</p>
<p>
In spite of its brightness, its presence is fleeting. The combined fact that we are actually seeing only reflected sunlight and that the ISS is in low-earth orbit means that there is only a brief window after sunset and before sunrise during which we can see this bright object in contrast against a dark sky. On top of that, the station’s orbit must pass over your geographic location within this window. This may seem to be a rare event when the rotation of the Earth and the steep inclination of the orbit are thrown into the mix. However, if you’re between the latitudes of -60&deg; and +60&deg;, then it happens and it happens predictably.
</p>
<p>
An overhead pass lasts roughly six minutes, with the station fading in as a steady point of light from one horizon, and fading out towards the opposite horizon. In this same time, the station tracks an arc of over 2,000 miles over North America. 
</p>
<p>
Bright as it is, an ISS pass itself can be wholly underwhelming. I’ve confirmed this in my on-campus attempts to photograph the station streaking over McGraw tower at Cornell University. My tripod-stabilized camera pointing skyward drew little more than confused looks from students passing by.
</p>
<img src="./images/ISSoverMcGraw.png" alt="Image of ISS flying high over McGraw Tower" style = "display: block; width:90%; margin-left: auto; margin-right:  auto; ">
<p>
To me, the sight of the ISS passing overhead affects a profound sense of humility and hope. Despite the retirement of the shuttle and the United States’ manned space program, the ISS has provided humankind with an uninterrupted 13 years in space. And in that passing point of light I see a future where humankind will look back at Earth and see that <a href="http://science.nasa.gov/science-news/science-at-nasa/2013/23jul_palebluedot/">pale blue dot</a> called home. A single speck of reflected sunlight harboring all that is and ever was of humanity and the only known life in the universe -- as the late Carl Sagan  famously described Earth in his commencement speech to the Cornell class of ’94. 
</p>
<p>
So look up every once in a while; next time you see a speck of light sailing through the night sky, take a pause from the stress and uncertainty of life. You are witness to one of humankind’s greatest accomplishments: an orbiting testament to our ingenuity, our progress, and our potential as a species.
</p>
<p style="text-align:center">Bryan Rolfe 2013 - bar72 [at] cornell.edu </p>
<p style="text-align:center">Two-line element data from: <a href="http://www.celestrak.com/">Celestrak</a> </p>

</div>


<script type="text/javascript" src="js/mercator.js"></script>

<script type="text/javascript">

	init();
    satellite();

    var colorSat     = 'rgba(0, 60, 255, 1)'; 
    var colorHome	 = 'rgba(0, 255, 60, 1)';
    var colorHorizon = 'rgba(50, 50, 50, 1)';
    var mapCalibration = false;

    function xCompare(a,b) {
		if (a.x < b.x)
			return -1;
		if (a.x > b.x)
			return  1;
		return 0;
	}

	function toRads(angle) {
		var rads = angle * (Math.PI) / 180;
		return rads;
	}

	function toXPixels(rads, center) {
		var pixels = center + (rads/Math.PI)*center*0.955 - 2;
		return pixels;
	}

	function toYPixels(rads, center, S) {
		var pixels = center - S*(rads)*center;
		return pixels;
	}

	function mercate(obj) {
		this.lon = obj.lon;
		this.lat = Math.log(Math.tan(Math.PI/4 + obj.lat/2));
	}

	function mercate2(lat) {
		return Math.log(Math.tan(Math.PI/4 + lat/2));
	}

	function compliment(theta) {
		theta = theta % (2*Math.PI);
		var sign = theta / Math.abs(theta);
		var theta2 = theta - sign * 2 * Math.PI; 
		var min = Math.min(Math.abs(theta), Math.abs(theta2));
		if (min == Math.abs(theta)) {
			return theta;
		} else {
			return theta2;
		}
	}

	function createFilledCircle(svg, x, y, r, color) {
		var circle = document.createElementNS(svgNS,'circle');
    	circle.setAttribute('cx', x);
    	circle.setAttribute('cy', y);
    	circle.setAttribute('r',r);
    	circle.setAttribute('fill',color);
    	circle.setAttribute('stroke','none');
    	circle.setAttribute('stroke-width','0px');

    	svg.appendChild(circle);
	}

	function createCircle(svg, x, y, r, color) {
		var circle = document.createElementNS(svgNS,'circle');
    	circle.setAttribute('cx', x);
    	circle.setAttribute('cy', y);
    	circle.setAttribute('r',r);
    	circle.setAttribute('fill','none');
    	circle.setAttribute('stroke',color);
    	circle.setAttribute('stroke-width','1px');

    	svg.appendChild(circle);
	}

	function createFilledPath(svg, points, color) {
		pathString = '';
		for (i = 0; i < points.length; i++) {
			if ( i==0 ) {
				pathString = pathString + 'M' + points[i].x + ',' + points[i].y;
			} else {
				pathString = pathString + " L" + points[i].x + ',' + points[i].y;
			}
		}

		//pathString = pathString + ' Z';

		var path = document.createElementNS(svgNS,'path');
		path.setAttribute('d', pathString);
		path.setAttribute('fill',color);
		path.setAttribute('fill-opacity',0.5);
		path.setAttribute('stroke', 'none');
		path.setAttribute('stroke-width', '0px');

		svg.appendChild(path);
	}

	function createPath(svg, points, color) {
		pathString = '';
		for (i = 0; i < points.length; i++) {
			if ( i==0 ) {
				pathString = pathString + 'M' + points[i].x + ',' + points[i].y;
			} else {
				pathString = pathString + " L" + points[i].x + ',' + points[i].y;
			}
		}

		//pathString = pathString + ' Z';

		var path = document.createElementNS(svgNS,'path');
		path.setAttribute('d', pathString);
		path.setAttribute('fill', 'none');
		path.setAttribute('stroke', color);
		path.setAttribute('stroke-width', '1px');

		svg.appendChild(path);
	}

	function checkPath(points, threshold) {
		var last = 0;
		for (var i = 0; i < ( points.length - 1 ) ; i++) {
			if ( Math.abs(points[i+1].x - points[i].x) > threshold ) {
				createPath(svg, points.slice(last,i), colorHorizon);
				last = i+1;
			}
		}
		createPath(svg, points.slice(last,i+1), colorHorizon);

		if ( Math.abs(points[i].x - points[0].x) <= threshold ) {
			createPath(svg, [points[i], points[0]], colorHorizon);
		}
	}

	function ready() {
		var container = document.getElementById("earthContainer");
		container.width = mapWidth;

		var svg   = document.getElementById("earth");
		var svgNS = svg.namespaceURI;

		svg.setAttribute("width", mapWidth);
		svg.setAttribute("height", mapHeight);

		//console.log(svg)
		var svgImage = document.createElementNS(svgNS, 'image');
		svgImage.setAttributeNS('http://www.w3.org/1999/xlink', 'href', "images/mercatorEarth.svg");
		svgImage.setAttribute('width', mapWidth);
		svgImage.setAttribute('height', mapHeight);
		svg.appendChild(svgImage);
	}

	function createCrosshair(x, y, color) {
		points = Array();
		var width = 10;
		var height = 10;

		points[0] = {x: x - width/2, y: y};
		points[1] = {x: x + width/2, y: y};
		createPath(svg, points, color);
		points[0] = {x: x, y: y - height/2};
		points[1] = {x: x, y: y + height/2};
		createPath(svg, points, color);
	}

	function createLatLonLines(color) {
		var latRange = [0, 15, 30, 45, 60, 75, -15, -30, -45, -60, -75];
		var lonRange = [0, 30, 60, 90, 120, 150, 180, -30, -60, -90, -120, -150];

		var points = Array();
		for (var i = 0; i < latRange.length; i ++ ) {
			points[0] = {x: 0, y: toYPixels(mercate2(toRads(latRange[i])), yZero, S)};
			points[1] = {x: mapWidth, y: toYPixels(mercate2(toRads(latRange[i])), yZero, S)};

			createPath(svg, points, color);
		}

		for (var i = 0; i < lonRange.length; i ++ ) {
			points[0] = {x: toXPixels(toRads(lonRange[i]), xZero), y: 0};
			points[1] = {x: toXPixels(toRads(lonRange[i]), xZero), y: mapHeight};

			createPath(svg, points, color);
		}
	}

	function clearSVGelements(svgDoc, name) {
		var temp = svgDoc.getElementsByTagName(name);

		for (i=temp.length-1; i >= 0; i--) {
			svgDoc.removeChild(temp[i]);
		}

	}

	// Function that returns a longitude for an object based on the current time
	//  we're converting from celestial coordinates to earth coordinates
	function correctLongitude(time, longitude) {
		var rotEquinox 	  = 2 * Math.PI * (time - vernalEquinox) / ( siderealFrac * dayMs);
    	var lonCorrection = equinoxLongi - rotEquinox % (2 * Math.PI) - Math.PI/2; 
    	// Although the Earth rotates eastward, the vernal equinox travels westward, since it is fixed in space relative to Earth,
    	//	this explains why there is a negative sign in front of the rotEquinox term (negative values are westward).
    	// The -pi/2 correction was added based on a hunch, and seems to produce a correct result, albeit I have only vague
    	//	ideas of why. I need to look into how the rotations of the orbit work, and I'm fairly positive I'll find a pi/2
    	//  correction elsewhere, that was used in the 3D projection.

    	if (lonCorrection >= 0) { 
    		var lonCorrected = (longitude + lonCorrection + Math.PI) % (2 * Math.PI) -  Math.PI;
    	} else {
    		var lonCorrected = (longitude + lonCorrection - Math.PI) % (- 2 * Math.PI) +  Math.PI;
    	}

	    return lonCorrected;	
	};

	function mercator() {

		// For orbits plotted over a mercator projection, the path of the orbit will be distorted from its 
	    // instanteous projection, since the earth is constantly rotating, and satellites don't orbit at the
	    // speed of infinity. 

		clearSVGelements(svg, "circle");
	    clearSVGelements(svg, "path");	

		if ( colorStyle == 'light' ) {
			 colorSat     = 'rgba(0, 60, 255, 1)'; 
    		 colorHorizon = 'rgba(50, 50, 50, 1)';
    		 colorLatLon  = 'rgba(0, 0, 0, 0.1)';
    		 colorInset   = 'rgba(0, 0, 0, 0.3)';
    		 colorHome	  = 'rgba(0, 190, 60, 1)';
		} else {
			 colorSat     = 'rgba(255, 0, 0, 1)'; 
    		 colorHorizon = 'rgba(255, 255, 255, 1)';
    		 colorLatLon  = 'rgba(255, 255, 255, 0.1)';
    		 colorInset   = 'rgba(255, 255, 255, 0.3)';
    		 colorHome	  = 'rgba(0, 190, 60, 1)';
		}

		canvasHolder.style.borderColor = colorInset;


	    //console.log("Active satellites:", activeSats.length, activeSats)
	    if (activeSats.length > 1) {
	    	//console.log("Sat Info: ", satInfo[activeSats[activeSats.length - 1]]);
	    	testOrbit = propagateOrbit(satInfo[activeSats[activeSats.length - 1]], 2*200, 2);
	    } else {
	    	testOrbit = propagateOrbit(satInfo[activeSats[0]], 2*200, 2);
	    }

	    var sunPixels = Array();
	    var sunMercator = Array();

	    for (i = 0; i < sunShadowLatLon.length; i ++ ) {
	    	sunPixels[i] = {x: toXPixels(sunShadowLatLon[i].lon, xZero),
	    					y: toYPixels(mercate2(sunShadowLatLon[i].lat), yZero, S)};
	    }

	    // After one sidereal day, the earth has rotated in such a way that the stars appear in the same place in the sky. 
	    // In other words, it is one complete rotation of the earth relative to the fixed background of the stars
	    // A solar day, 24h, is one rotation such that the Sun appears in the same position in the sky (disregarding elevation above the horizon)
	    // The important implication is that, after any integer number of sidereal days, the earth-longitude where the equinox occured
	    //	will exactly line up with the celestial vernal equinox direction. 

	    rotEquinox = tNowMs - vernalEquinox;
		rotEquinox = 2 * Math.PI * rotEquinox/( siderealFrac * dayMs);

		//                   |------------------- Vernal Equinox Earth-Longitude
		//                   |              |---- Sidereal rotations (radians) from vernal equinox
		//                   |              |		if this value is an integer, than the modulo = 0, and the longitude is the equinox longitude
	    lonCorrection = equinoxLongi - rotEquinox % (2 * Math.PI);
	    lonCorrection = compliment(lonCorrection);

	    var vernalPath = Array();
	    vernalPath[0] = {x: toXPixels( lonCorrection, xZero), y: 0};
	    vernalPath[1] = {x: toXPixels( lonCorrection, xZero), y: mapHeight};
	    createPath(svg, vernalPath, 'rgba(255, 0, 0, 0.4)');
		
	    var lonCorrected = 0;

	    //console.log("Long. Correction:", 180 * (lonCorrection % ( 2*Math.PI)) / Math.PI)


	    // Horizon points for current satellite position

	    horizonPixels = Array()
	    for (i = 0; i < horizonLatLon.length; i ++ ) {
	    	lonCorrected = correctLongitude(testOrbit[0].t, horizonLatLon[i].lon);
	    	horizonPixels[i] = {x: toXPixels(lonCorrected, xZero),
	    					y: toYPixels(mercate2(horizonLatLon[i].lat), yZero, S)};

	    	//createFilledCircle(svg, toXPixels(lonCorrected, xZero), toYPixels(mercate2(horizonLatLon[i].lat), yZero, S), 1, 'grey');
	    }

	    // Satellite path projection and longitude corrections
	    orbitPixels = Array();
	    latLon = Array();

	    for (i = 0; i < testOrbit.length; i ++ ) {

	    	lonCorrected = correctLongitude(testOrbit[i].t, Math.atan2(testOrbit[i].x, testOrbit[i].z));
	    	//console.log("Longitude correction: ", rotEquinox, rotEquinox % (2 * Math.PI), lonCorrection, lonCorrection - Math.PI, (lonCorrection - Math.PI) % (-2*Math.PI), (lonCorrection - Math.PI) % (-2*Math.PI) + Math.PI)
	    	//console.log("Orbit time:", testOrbit[i].t, vernalEquinox)

	    	d = norm([testOrbit[i].x, testOrbit[i].y, testOrbit[i].z]);
	    	latLon[i] = {lon: lonCorrected ,
						 lat: Math.PI/2 - Math.acos(testOrbit[i].y / d )};

			orbitPixels[i] = {x: toXPixels(latLon[i].lon, xZero),
	    					y: toYPixels(mercate2(latLon[i].lat), yZero, S)};

	    	createFilledCircle(svg, orbitPixels[i].x, orbitPixels[i].y, 1, colorSat);

	    }

	    // Calibration Points
	    if (mapCalibration) {
		    createCircle(svg, toXPixels(toRads(-155.674), xZero), toYPixels(mercate2(toRads(19.67)), yZero, S), 3, 'red');
		    createCircle(svg, toXPixels(toRads(130.609), xZero), toYPixels(mercate2(toRads(31.388)), yZero, S), 3, 'red');
		    createCircle(svg, toXPixels(toRads(142.298), xZero), toYPixels(mercate2(toRads(-11.009)), yZero, S), 3, 'red');
		    createCircle(svg, toXPixels(toRads(-65.563), xZero), toYPixels(mercate2(toRads(-54.714)), yZero, S), 3, 'red');
		    createCircle(svg, toXPixels(toRads(19.5325), xZero), toYPixels(mercate2(toRads(-34.561)), yZero, S), 3, 'red');
		    createCircle(svg, toXPixels(toRads(-7.0224), xZero), toYPixels(mercate2(toRads(62.193)), yZero, S), 3, 'red');
		    createCircle(svg, toXPixels(toRads(-36.512), xZero), toYPixels(mercate2(toRads(-54.408)), yZero, S), 3, 'red');
		    createCircle(svg, toXPixels(toRads(146.477), xZero), toYPixels(mercate2(toRads(-43.456)), yZero, S), 3, 'red');
		}

		createLatLonLines(colorLatLon);

	    // Sort by x-position (pixels) -- sorting by longitude, essentually.
	    sunPixels.sort(xCompare);
	    if ( colorStyle == 'dark' ) {
	    	createPath(svg, sunPixels, 'white');
	    }
	    sunPixels[sunPixels.length] = {x: toXPixels(Math.PI,xZero), y: toYPixels(mercate2(toRads(-76)), yZero, S)};
	    sunPixels[sunPixels.length] = {x: toXPixels(-Math.PI,xZero), y: toYPixels(mercate2(toRads(-76)), yZero, S)};
	    sunPixels[sunPixels.length] = {x: sunPixels[0].x, y: sunPixels[0].y};
	    createFilledPath(svg, sunPixels, 'black');

	    createCrosshair(orbitPixels[0].x, orbitPixels[0].y, colorSat);

	    //horizonPixels.sort(xCompare);
	    //createPath(svg, horizonPixels);
	    //console.log(horizonPixels);
	    checkPath(horizonPixels, xZero);

	    createCrosshair(toXPixels(toRads(homeLongi), xZero), toYPixels(mercate2(toRads(homeLat)), yZero, S), colorHome);

	    
	}

	function update() {
		//init();
    	satellite();
		mercator();
		timer = setTimeout(update, tUpdate);
	}

    var svg   = document.getElementById("earth");

    var svgNS = svg.namespaceURI;
    var mapWidth  = 800;
    var mapHeight = (443/600) * mapWidth;
    var xZero    = mapWidth/2;
    var offset = (38.13/443) * mapHeight;
    var yZero = mapHeight/2 + offset;
    var S = (1/2.88247); // Scaling factor for the mercator projection based on the specific SVG image I used
    
    ready();
    mercator();

    //console.log("Tangent testing: ", Math.atan2(0,1), Math.atan2(1,0), Math.atan2(0, -1), Math.atan2(-1,0))

    // TEMPORARY
    if  ( updateSw ) {
	    var timer = setTimeout(update, tUpdate) 
	}
	// TEMPORARY

    // Event Listeners
	// Need to come back and edit this section for browser compatability, especially older browers like IE5-8.
    
    var closeInset   = document.getElementById("close");
    var canvasHolder = document.getElementById("canvasHolder");

    closeInset.onclick = function() {
    	draw3D = ! draw3D;
    	console.log(draw3D);
    	if (draw3D == true) {
    		canvasHolder.style.display = 'block';
    		console.log("visible")
    		this.text = "X";
    		init();
    		satellite();
    	} else {
    		canvasHolder.style.display = 'none';
    		console.log("not visible")
    		this.text = "O";
    		init();
    		satellite();
    	}
    }
	dropdownMenu.addEventListener('change', function(e) { 
			if (e.target.value >= 0 && activeSats.find(e.target.value) == -1) {
				activeSats.remove(currentSat);
				currentSat = parseInt(e.target.value); 
				activeSats.remove(currentSat); // If it exists, remove it to avoid duplicates
				activeSats.push(currentSat); 
				satInfo[currentSat] = new setSatellite(TLEobjs[currentSat]); 

				//a = document.createElement('a');
				//var linkText = document.createTextNode(satInfo[currentSat].objectName);
				//a.appendChild(linkText);
				//a.href = "#";
				//a.id = currentSat;
				//a.onclick = function () { 
				//	activeSats.remove(parseInt(this.id));
				//	currentSat = activeSats[activeSats.length - 1];
				//	init(); 
				//	satellite();
				//	mercator();

				//	this.parentNode.removeChild( this ) 
				//};

				//elementActiveSats.appendChild(a);

				init(); 
				satellite();
				mercator();
			}
		}, false);

	elementULongi.addEventListener('change', function(e) {init()}, false);
	elementULat.addEventListener('change', function(e) {init()}, false);
	elementTimeForward.addEventListener('change', function(e) {init()}, false);



</script>
</html>