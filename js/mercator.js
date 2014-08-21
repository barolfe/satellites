//var elementSatRA  = document.getElementById("satRA");
//var elementSatDEC = document.getElementById("satDEC");
var elementTime   	 	= document.getElementById("currentTime");
var elementEpoch  	 	= document.getElementById("epochTime");
//var elementName   	 	= document.getElementById("objectName");
var elementEccentric 	= document.getElementById("eccentricAnomoly");
var elementInclination 	= document.getElementById("inclination");
var elementaoPerigee	= document.getElementById("aoPerigee");
var elementraAscending	= document.getElementById("raAscending");
var elementMeanAnomoly	= document.getElementById("meanAnomoly");
var elementMeanMotion	= document.getElementById("meanMotion");
var elementEccentricity = document.getElementById("eccentricity");
var elementVelocity		= document.getElementById("velocity");
var elementAltitude		= document.getElementById("altitude");
var elementISSspeed		= document.getElementById("ISSspeed");
var elementISSalt		= document.getElementById("ISSalt");
//var elementLongiLat		= document.getElementById("LongiLat");
var elementULongi 		= document.getElementById("uLongi");
var elementULat			= document.getElementById("uLat");
var elementTimeForward  = document.getElementById("timeForward");
var elementActiveSats   = document.getElementById("activeSatellites");
//var tUpdate = 1000; // ms, frequency to update satellite calculations and drawing
var updateSw = true;
var tUpdate = 1000; // ms, frequency to update satellite calculations and drawing
var draw3D = true;
var siderealFrac = 0.99726958;
var dayMs = 86400000;
var timeForward = 0.0;

function getUsersCoordinates() {
	if (navigator.geolocation) {
		console.log('success');
		navigator.geolocation.getCurrentPosition(assignCoords);
		return true; 
	} else {
		console.log("geolocation is not supported...");
		return false;
	}
}

function assignCoords(position) {
	homeLat   = position.coords.latitude;
	homeLongi = position.coords.longitude;
	elementULat.value = homeLat.toFixed(4);
	elementULongi.value = homeLongi.toFixed(4);
	init();
	console.log(position.coords);
}


function toggle_visibility(id) 
{   var e = document.getElementById(id);
   if(e.style.display == 'block')
      e.style.display = 'none';
   else
      e.style.display = 'block';
}

function toggle_color() {
   var e = document.body;
   if(e.style.background == 'rgb(255, 255, 255)') {
      e.style.background = 'rgb(0, 0, 0)';
      e.style.color 	 = 'rgb(200, 200, 200)';
      colorStyle 		 = 'dark';
   } else {
      e.style.background = 'rgb(255, 255, 255)';
      e.style.color 	 = 'rgb(40, 40, 40)';
      colorStyle 		 = 'light';
   }
}

//Array prototype function to remove a specified element(s) from an array
Array.prototype.remove = function() {
	var what, a = arguments, L = a.length, ax;
	while (L && this.length) {
		what = a[--L];
		while  ((ax = this.indexOf(what)) !== -1) {
			this.splice(ax, 1);
		}
	}
	return this;
}

Array.prototype.find = function(a) {
	var index = -1;
	for (i=0; i < this.length; i++) {
		if (a == this[i]) {
			index = i;
		}
	}

	return index;
}

function cross(u, v) {
	var U = Array();

	U[0] = u[1]*v[2] - u[2]*v[1];
	U[1] = u[2]*v[0] - u[0]*v[2];
	U[2] = u[0]*v[1] - u[1]*v[0];

	return U;
}

function dot(u, v) {
	
	var d = u[0]*v[0] + u[1]*v[1] + u[2]*v[2];

	return d;
}

function norm(u) {
	return Math.sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
}

function createVector(pStart, pEnd, numPts, radius, rgb) {
	var lineObj = Array();
	var Dx = pEnd[0] - pStart[0]; var Dy = pEnd[1] - pStart[1]; var Dz = pEnd[2] - pStart[2];
	var x0 = pStart[0]; var y0 = pStart[1]; var z0 = pStart[2];

	for (var i = 0; i < numPts; i++) {
		lineObj[i] = {x: (i/numPts) * Dx + x0,  y: (i/numPts) * Dy + y0, z: (i/numPts) * Dz + z0, r: radius, red: rgb[0], green: rgb[1], blue: rgb[2]};
	}

	this.pts = lineObj;
	this.D = [Dx, Dy, Dz];
	this.v0 = [x0, y0, z0];
}

// Return latitude and logitude information for an object array of x, y, z data in units pixels
function returnLatLon(xyzObj) {
	var latLonObj = Array();
	for (var i = 0; i < xyzObj.length; i++) {
	latLonObj[i] = {lon: Math.atan2(xyzObj[i].x, xyzObj[i].z),
						  		  lat: Math.PI/2 - Math.acos(xyzObj[i].y / ( (earthRadius*1000) / metersToPix) )};
	}
	return latLonObj;
}

function getEccentricAnomoly(satInfo) {
	// This has to be done numerically
	var Enew = satInfo.totalAnomoly; // Initial Guess
	var diffThreshold = 2 * Math.PI/3600; // accurate to within a tenth of a degree... for now
	var diff = 10*diffThreshold;

	while (diff > diffThreshold) {
		var Eold = Enew;
		var zero = Enew - satInfo.eccentricity*Math.sin(Enew) - satInfo.totalAnomoly;
		var slope = 1 - satInfo.eccentricity*Math.cos(Enew);
		Enew = Eold - zero/slope; 						// Enew = Eold - f(Eold)/f'(Eold)
		diff = Math.abs(Enew - Eold);
	}
	
	return Enew; 
}


// createOrbit will return an array of points that outline the path of the orbit fitting the 
// 	provided parameters. This function returns points in the celestial coordinate system. 
//	This function does NOT include the position of the satellite at any given time, only
//	the orbit in space. 
function createOrbit(eccentricity, aoPerigee, inclination, raAscending, radius, numCircle) {
	// Orbit points
	var orbit 		= Array();
	var	aEllipse 	= 1;
	var bEllipse 	= aEllipse * Math.sqrt(1 - eccentricity*eccentricity);
	var f 			= Math.sqrt(aEllipse*aEllipse - bEllipse*bEllipse);  // distance from center of ellipse to focci

	aEllipse = radius * aEllipse;
	bEllipse = radius * bEllipse;
	f        = radius * f;

	for (i=0; i < numCircle; i++) {

		theta = (i/numCircle)*Math.PI*2;
		orbit[i] = {x: aEllipse*Math.cos(theta) + f, y: 0, z: bEllipse*Math.sin(theta), r: 1, red: rgb[0], green: rgb[1], blue: rgb[2]};

		rotateParticle(orbit[i], 0, aoPerigee - 3 * Math.PI/2, 0)
		rotateParticle(orbit[i], 0, 0, inclination);
		rotateParticle(orbit[i], 0, raAscending, 0);

	}

	return orbit;
}

// The purpose of this function will be to provide orbit data for the purposes of an earth-latitude/longitude projection
//	This function will behave similarly to the createOrbit function, with the key different being that each point provided
//	will have a timestamp associated with it that will/can be used to correctly project the satellite onto Earth.

function propagateOrbit(uSatInfo, numCircle, N) {

	var orbitalPoints = Array();

	var eccentricity = uSatInfo.eccentricity;
	var aoPerigee    = uSatInfo.aoPerigee;
	var inclination  = uSatInfo.inclination;
	var raAscending  = uSatInfo.raAscending;
	var epoch        = uSatInfo.epoch;
	var meanAnomoly  = uSatInfo.meanAnomoly;
	var meanMotion   = uSatInfo.meanMotion;
	var aEllipse     = uSatInfo.aEllipse;
	var bEllipse     = uSatInfo.bEllipse;
	var f            = uSatInfo.f;
	var epochObj	 = uSatInfo.epochObj;

	var rScaleFactor = uSatInfo.a / satInfo[currentSat].a;

	var tnow = new Date();
	var tNowMs = tnow.getTime() + timeForward*60*60*1000;


	aEllipse = aEllipse * rScaleFactor;
 	bEllipse = bEllipse * rScaleFactor;
 	f        = f * rScaleFactor;

 	t =  tNowMs  - epoch; 		 // time elapsed since epoch, milliseconds
	t = t/(24 * 60 * 60 * 1000); // time elapsed since epoch, days

 	var satObjTemp = {eccentricity: eccentricity, totalAnomoly: meanAnomoly + t * meanMotion};

 	var timeEnd = N * 2 * Math.PI / meanMotion; // time (days) required N orbits
 	var timeEnd = timeEnd * 24 * 60 * 60 * 1000; // time (milliseconds)

 	//console.log("Time to orbit: ", timeEnd / (1000 * 60), " minutes");

 	var delTime = timeEnd / numCircle; // time step interval such that numCircle * delT = timEnd

 	var tMs = 0;

 	var tFromEpoch =  tNowMs  - epoch; 		 // time elapsed since epoch, milliseconds
 	//console.log("Orbit start time:", tnow);

	for (var i = 0; i < numCircle; i++) {

		t = tFromEpoch + i * delTime;		 // time propagated forward, from epoch, milliseconds
		tMs = t;
		t = t / (24 * 60 * 60 * 1000);       // time propagated forward, from epoch, days
		tMs = tMs - tFromEpoch;				 // time propagated forward, starting from 0 = current time
		tMs = tMs + tNowMs;					 // absolute time, propagated forward, UNIX time.

		satObjTemp.totalAnomoly = meanAnomoly + t * meanMotion;
		E = getEccentricAnomoly(satObjTemp);

		orbitalPoints[i] = {x: aEllipse*Math.cos(Math.PI - E) + f,
						 	y: 0,
						 	z: bEllipse*Math.sin(Math.PI - E),
						 	t: tMs,
						 	r: 1, 
						 	red: 0,
						 	green: 0,
						 	blue: 255};

		rotateParticle(orbitalPoints[i], 0, aoPerigee +  Math.PI/2, 0)
		rotateParticle(orbitalPoints[i], 0, 0, inclination);
		rotateParticle(orbitalPoints[i], 0, raAscending, 0);
	}

	return orbitalPoints;

}

// This function returns an array of lat/lon objects marking the visible horizon of a given satellite at a given time
function visibleHorizon(satObj, alt) {
	// A satellite's visible horizon is calculated from simple geometry. 
	//   The horizon is defined by line-of-sight, meaning lines connecting the satellite through tangets of the Earth
	//   mark this horizon. Connecting tangets, to the center of the Earth (radius R), and these two points to the 
	//   satellite, form a right angle triangle. Theta, the angle between center-tangent and center-sat, is calculated by:
	//   	theta = acos(R / ( R + alt ));
	//   Then the unknown leg of the triangle, tangent-sat (c), can be calcuated:
	//      c = ( R+ alt )*sin(theta);
	//   Finally, we form a new right triangle, formed by tanget-sat, sat-point, point-tanget -- where point is the intersection
	//   of the line R + h (center-sat) and a perpendicular line that intersects the tanget point. The hypotonuse of this triangle
	//   is already known, as well as the angle theta, (similar triangles), so we can readily calculate the length of point-tangent (a)
	//      a = c * cos(theta)
	//   Thus:
	//   	a = ( R + alt ) * sin(acos( ( R / ( R + alt ) ) * cos(acos( R / ( R + alt ) ))
	//      a = ( R + alt ) * sqrt( 1 - ( R / ( R + alt )) ^ 2 ) * R / ( R + alt )
	//

	var theta = Math.acos(earthRadius * 1000 / ( earthRadius*1000 + alt ));
	var c     = (earthRadius * 1000 + alt) * Math.sin( theta );
	var a 	  = c * Math.cos(theta); // Distnace from perpendicular point to tangent
	var b 	  = c * Math.sin(theta); // Distance from perpendicular point to satellite

	var satVect       = [satObj.x * metersToPix, satObj.y * metersToPix, satObj.z * metersToPix]; // meters
	var satVectLength = norm(satVect);
	var perpVect      = Array();

	for (var i = 0; i < 3; i++) { perpVect[i] = satVect[i] * (satVectLength - b)/(satVectLength); }

	var U = Array();
	var V = Array();

	for (var i = 0; i < 3; i++) { U[i] = perpVect[i]; } // normalize a vector pointing to the satellite
    V = [U[0], 0, 0];
	var W = cross(U, V);
	V = cross(U, W);

	// V, W, are now two vectors normal to the vector pointing to the satellite. 

	var VLength = norm(V);
	var WLength = norm(W);
	for (var i = 0; i < 3; i++) { 
		V[i] = (V[i] / VLength) * a;
		W[i] = (W[i] / WLength) * a;
	} 	

	// V, W, are now normalized to the length of a, the radius of the horizon circle

	var horizonPoints = Array();
	for (var i = 0; i < 100; i++) {
		angle = (i/100) * 2 * Math.PI;

		horizonPoints[i] = {x: (V[0]*Math.cos(angle) + W[0]*Math.sin(angle) + perpVect[0]) / metersToPix, 
						y: (V[1]*Math.cos(angle) + W[1]*Math.sin(angle) + perpVect[1]) / metersToPix,
						z: (V[2]*Math.cos(angle) + W[2]*Math.sin(angle) + perpVect[2]) / metersToPix,
						r: 1, red: 200, green: 0, blue: 0};
	}

	return horizonPoints;

}


function setSatellite(satObj) {
	// All measurements are in the Equitorial Coordinate System
	// 	The epoch stuff is a bit confusing. The current date is in reference to midnight on December 31st, not January 1st. This is correct, but why?
	this.epoch = Date.UTC(satObj.year, 0, 0, 0, 0, 0, 0) + satObj.epoch * 24 * 60 * 60 * 1000; 
	this.epochObj = new Date(this.epoch);

	this.objectName  = satObj.name;
	this.inclination = 2*Math.PI*(satObj.inclination /360);				// radians
	this.raAscending = 2*Math.PI*(satObj.raAscending /360) + Math.PI/2;	// radians, the pi/2 is necessary because of where I chose the ascending node to be on the ellipse
	this.eccentricity = satObj.eccentricity;								// no units
	this.aoPerigee   = 2*Math.PI*(satObj.aoPerigee   /360);				// radians

	this.meanAnomoly = 2*Math.PI*(satObj.meanAnomoly /360);				// radians (not technically an angle though)
	this.meanMotion  = 2*Math.PI * satObj.meanMotion;					// rads/day

	// Elliptical orbits (all real orbits are somewhat elliptical)
	//	e - eccentricity
	//		e = sqrt(1 - (b/a)^2); e^2 = 1 - (b/a)^2; (b/a) = sqrt(1 - e^2)
	//	M = E - e * sin(E), E is the eccentric anomoly;  no closed form solution!
	//	ellipse: x^2/a^2 + y^2/b^2 = 1
	//	If we know E, a, and b:
	//		x = a * cos(E);  y = b * sin(E)
	//	Focci -- Earth is at one of the two focci 
	//		f = sqrt( a^2 - b^2 ) = distance of focci from center of ellipse

	this.a 			= Math.pow( ( G * earthMass ) / Math.pow(this.meanMotion / (86400) , 2 ), 1.0/3.0 ); 	// m

	this.aEllipse 	= 1;
	this.bEllipse 	= this.aEllipse * Math.sqrt(1 - this.eccentricity*this.eccentricity);
	this.f 			= Math.sqrt(this.aEllipse*this.aEllipse - this.bEllipse*this.bEllipse);  			// distance from center of ellipse to focci

	this.earthR 		= refRadius * (1000 * earthRadius / this.a);
	
	this.aEllipse 	= refRadius * this.aEllipse;
	this.bEllipse 	= refRadius * this.bEllipse;
	this.f        	= refRadius * this.f;

	this.maxDepth = this.aEllipse + this.f;

	//console.log(this)
}

function doOnce() {
}

function init() {

	// Set color style for the page
	if (colorStyle == 'dark') {
		rgb = [255, 255, 255];
		canvasBG = 'rgb(0, 0, 0)';
	} else {
		rgb = [10 , 10 , 10 ];
		canvasBG = 'rgb(255, 255, 255)';
	}

	// Get time to offset the calculations by
	timeForward = elementTimeForward.value;

	// Get user's coordinates
	homeLat		= elementULat.value;
	homeLongi	= elementULongi.value;

	// Assign variables for convenience -- This is a temporary fix
	var eccentricity = satInfo[currentSat].eccentricity;
	var aoPerigee    = satInfo[currentSat].aoPerigee;
	var inclination  = satInfo[currentSat].inclination;
	var raAscending  = satInfo[currentSat].raAscending;
	var epoch        = satInfo[currentSat].epoch;
	var meanAnomoly  = satInfo[currentSat].meanAnomoly;
	var meanMotion   = satInfo[currentSat].meanMotion;
	var aEllipse     = satInfo[currentSat].aEllipse;
	var bEllipse     = satInfo[currentSat].bEllipse;
	var f            = satInfo[currentSat].f;

	var earthR       = satInfo[currentSat].earthR;

	metersToPix = satInfo[currentSat].a /refRadius; 

	// Orbit points
	var rScaleFactor = satInfo[currentSat].a / satInfo[currentSat].a;
 	satOrbit[currentSat] = createOrbit(eccentricity, aoPerigee, inclination, raAscending, rScaleFactor * refRadius, numCircle);

 	// Make sure the satellite orbit is colored according to the current color style -- Come back to this

 	aEllipse = aEllipse * rScaleFactor;
 	bEllipse = bEllipse * rScaleFactor;

	earthObj = {x: 0, y: 0, z: 0, r: earthR, red: 0, green: 150, blue: 255};

	tTemp = new Date()
	tnow = tTemp.getTime();

	t =  tnow - epoch; 	// time elapsed since epoch, milliseconds

	t = t/(24 * 60 * 60 * 1000); // time elapsed since epoch, days


	satInfo[currentSat].totalAnomoly = meanAnomoly + t * meanMotion;

	// For ellipical orbits, we need the eccentric anomoly. 
	E = getEccentricAnomoly(satInfo[currentSat]);

	satObj[currentSat] = {x: aEllipse*Math.cos(Math.PI - E), y: 0, z: bEllipse*Math.sin(Math.PI - E), r: 4, red: 255, green: 0, blue: 0 };

	rotateParticle(satObj[currentSat], 0, aoPerigee + Math.PI/2, 0)
	rotateParticle(satObj[currentSat], 0, 0, inclination);
	rotateParticle(satObj[currentSat], 0, raAscending, 0);

	// For veloctiy calculation
	pold.x = satObj[currentSat].x; pold.y = satObj[currentSat].y; pold.z = satObj[currentSat].z;
	told = tnow; // ms

	for (i=0; i < numLongi; i++) {

		theta = (i/numLongi) * Math.PI * 2;
		longiCircle[i] = {x: earthR*Math.cos(theta), y:  earthR*Math.sin(theta), z: 0, r: 1, red: 255, green: 255, blue: 255};

	}

	var longiCnt = 0;
	for (i=0; i < 36; i++) {
		longitude = 2 * Math.PI * (i / 36);
		for (k=0; k < numLongi; k++) {
			longiObj[longiCnt] = {x: longiCircle[k].x, y: longiCircle[k].y, z: longiCircle[k].z, r: 1, red: 255, green: 255, blue: 255};
			rotateParticle(longiObj[longiCnt], 0, longitude, 0);
			longiCnt ++;
		}	

	}

	homeObj[0] = {x: earthR, y:0, z:0, r:4, red: 0, green: 200, blue: 0 };

	// Now we're working in equitorial coordinates. We'll use a series of rotations to positon our home particle.
	rotateParticle(homeObj[0], 0, 0, 2 * Math.PI * homeLat/360)						//Position latitude first so we can keep the calculations simple
	rotateParticle(homeObj[0], 0, (2 * Math.PI *homeLongi/360 - equinoxLongi), 0);		//Position longitude at vernal equinox.

	elementInclination.innerHTML 	= (inclination/(2*Math.PI)*360).toFixed(4) + "&deg;";
	elementaoPerigee.innerHTML	 	= (aoPerigee/(2*Math.PI)*360).toFixed(4)	+ "&deg;";
	elementraAscending.innerHTML 	= (raAscending/(2*Math.PI)*360).toFixed(4)	+ "&deg;";
	elementMeanAnomoly.innerHTML 	= (meanAnomoly/(2*Math.PI)*360).toFixed(4) + "&deg;";
	elementMeanMotion.innerHTML	 	= (meanMotion/(2*Math.PI)).toFixed(4)		+ " orbits/day";
	elementEccentricity.innerHTML 	= eccentricity.toFixed(6);

	// Now we're going to calculate the orientation of the Earth at the current date, relative to the vernal equinox.
	//	We have to be careful to use sidereal days here, and not solar days, as we want to leave the vernal equinox
	//	aligned to the positive x-axis for convenience. One sidereal day, on average, is 23 hours, 56 minutes and 4 
	//	seconds, or 0.99726958 solar days.
}

function zCompare(a,b) {
	if (a.z < b.z)
		return -1;
	if (a.z > b.z)
		return  1;
	return 0;
}

function renderObjs(context, objs, view) {
	// Render array of objs, but first create a copy to work with as we don't want to lose the ordering!
	//var tempObjs = objs.clone(); 					// clone() is an expensive routine
	
	// simpe view rotation on all particles
	for (i=0; i < objs.length; i++) {
		rotateParticle(objs[i], 0, view.y, 0)
		rotateParticle(objs[i], view.x, 0, 0)
	}

	// sort particles by z for proper depth rendering
	objs.sort(zCompare);

	for (i=0; i < objs.length; i++) {
		drawParticle(objs[i].x, objs[i].y, objs[i].z, objs[i].r, objs[i].red, objs[i].green, objs[i].blue);
	}

	// reutrn our object array to its original state
	objs.sort(function(a,b) { return (a.key - b.key) });

	for (i=0; i < objs.length; i++) {
	 	rotateParticle(objs[i], -view.x, 0, 0)
	 	rotateParticle(objs[i], 0, -view.y, 0)	
	}

}

function drawParticle(xp, yp, zp, r, red, green, blue) {

	xp = xp + xc;
	// the y-value is flipped because the canvas coordinate
	// system has its orgin in the top-left
	yp = -yp + yc; 

	alpha = (zp + maxDepth + 10)/(2*maxDepth + 10);

	context.fillStyle = "rgba(" + red + "," + green + "," + blue + "," + alpha +")";
	context.beginPath();
	context.arc(xp, yp, r, 0, 2*Math.PI);
	context.closePath();
	context.fill();

}

function rotateParticle(pObj, xtheta, ytheta, ztheta) {

	// function that rotates the x, y, z coordinates of the passed object
	// about an angle ztheta, then ytheta, and then xtheta
	// all rotations follow the right-hand rule

	var xOld = pObj.x;
	var yOld = pObj.y;
	var zOld = pObj.z;

	pObj.x = xOld*Math.cos(ztheta) - yOld*Math.sin(ztheta);
	pObj.y = xOld*Math.sin(ztheta) + yOld*Math.cos(ztheta);

	xOld = pObj.x;
	yOld = pObj.y;
	zOld = pObj.z;

	pObj.x =  xOld*Math.cos(ytheta) + zOld*Math.sin(ytheta);
	pObj.z = -1*xOld*Math.sin(ytheta) + zOld*Math.cos(ytheta);

	xOld = pObj.x;
	yOld = pObj.y;
	zOld = pObj.z;

	pObj.y = yOld*Math.cos(xtheta) - zOld*Math.sin(xtheta);
	pObj.z = yOld*Math.sin(xtheta) + zOld*Math.cos(xtheta);

}

// Satellite Position
function satellite() {

	// Get current date object to compare determine how many milliseconds have ellapsed in the day
	tnow = new Date();
	t0   = new Date(tnow.getUTCFullYear(), tnow.getUTCMonth(), tnow.getUTCDate(), 0, 0, 0, 0);

	tNowMs = tnow.getTime() + timeForward*60*60*1000;

	var earthR       = satInfo[currentSat].earthR;

	// Loop over the currently displayed satellites, and update their positions for the present time. 
	for (i = 0; i < activeSats.length; i++) {
		var eccentricity = satInfo[activeSats[i]].eccentricity;
		var aoPerigee    = satInfo[activeSats[i]].aoPerigee;
		var inclination  = satInfo[activeSats[i]].inclination;
		var raAscending  = satInfo[activeSats[i]].raAscending;
		var epoch        = satInfo[activeSats[i]].epoch;
		var meanAnomoly  = satInfo[activeSats[i]].meanAnomoly;
		var meanMotion   = satInfo[activeSats[i]].meanMotion;
		var aEllipse     = satInfo[activeSats[i]].aEllipse;
		var bEllipse     = satInfo[activeSats[i]].bEllipse;
		var f            = satInfo[activeSats[i]].f;
		var epochObj	 = satInfo[activeSats[i]].epochObj;

		var rScaleFactor = satInfo[activeSats[i]].a / satInfo[currentSat].a;
		var tempMaxDepth = rScaleFactor * satInfo[activeSats[i]].maxDepth;

		if (tempMaxDepth > maxDepth) maxDepth = tempMaxDepth;

		aEllipse = aEllipse * rScaleFactor;
	 	bEllipse = bEllipse * rScaleFactor;
	 	f        = f * rScaleFactor;

	 	t =  tNowMs  - epoch; 		 // time elapsed since epoch, milliseconds
		t = t/(24 * 60 * 60 * 1000); // time elapsed since epoch, days

	 	satInfo[activeSats[i]].totalAnomoly = meanAnomoly + t * meanMotion;

		// For ellipical orbits, we need the eccentric anomoly. 
		E = getEccentricAnomoly(satInfo[activeSats[i]]);

		// The negative signs have to do with how the orbit is defined in the TLE Data
		//	and how I am defining the satellite position on a circle
		satObj[activeSats[i]].x = aEllipse*Math.cos(Math.PI - E) + f;
		satObj[activeSats[i]].y = 0;
		satObj[activeSats[i]].z = bEllipse*Math.sin(Math.PI - E);

		rotateParticle(satObj[activeSats[i]], 0, aoPerigee + Math.PI/2, 0)
		rotateParticle(satObj[activeSats[i]], 0, 0, inclination);
		rotateParticle(satObj[activeSats[i]], 0, raAscending, 0);
	}

	var tday = tNowMs - t0.getTime();
	var dayFraction = tday/( siderealFrac * dayMs); // fraction of mean sidereal day, not solar day

	rotEquinox = tNowMs - vernalEquinox;
	rotEquinox = 2 * Math.PI * rotEquinox/( siderealFrac * dayMs);

	//////////////////////
	// SUN CALCULATIONS //
	//////////////////////

	// Must be very careful about confusing earth-based lattitude and logitudes with celestial lattitudes and longitudes
	// (don't panic)
	solarLapsed = tNowMs - vernalEquinox;
	//sunLongi 	= Math.PI * 2 * (249/360);

	// By definition, one day -- dayMs seconds -- is a single rotation that places the sun at the same longitude
	// Hence, an integer number of days will place the sun at the equinox longitude, anything else will be the 
	//  the fraction of the day multiplied into 360 degrees subtracted from the position of the sun at equinox.
	//  The negative sign makes sense because the term (tNowMs - vernalEquinox) is measured westward, because
	//  the earth rotates westward. 

	sunLat		= 2 * Math.PI * (23.4/360) * Math.sin(2 * Math.PI * solarLapsed/(dayMs * 365));
	// Sun's longitude over the Earth -- this used for the mercator projection
	sunLongi = 360 * equinoxLongi / (2 * Math.PI) - 360 * ( solarLapsed / dayMs) % 360;
	// Sun's longitude relative to vernal equinox (celestial longitude) -- this is used for the 3D representaiton
	sunLongiVernal = equinoxLongi + 2 * Math.PI * solarLapsed / (dayMs * 365) - Math.PI; 
	sunLongi = sunLongi * Math.PI / 180;


	xsun = earthR * Math.cos(sunLat) * Math.sin(sunLongi);
	zsun = earthR * Math.cos(sunLat) * Math.cos(sunLongi);
	ysun = earthR * Math.sin(sunLat);

	xsunV = earthR * Math.cos(sunLat) * Math.sin(sunLongiVernal);
	zsunV = earthR * Math.cos(sunLat) * Math.cos(sunLongiVernal);
	ysunV = earthR * Math.sin(sunLat);

	// Normalize the vectors that point to the sun 
	sunLen = Math.sqrt(xsun * xsun + zsun * zsun + ysun * ysun);
	xsunN = xsun / sunLen;
	zsunN = zsun / sunLen;
	ysunN = ysun / sunLen;
	
	sunV = [xsunN, ysunN, zsunN];
	refV = [xsunN, 0, 0];
	greatCircleU = Array();

	// Compute cross product of sun vector and an arbitrary reference vector
	greatCircleU = cross(sunV, refV);
	greatCircleV = cross(sunV, greatCircleU);
	normU = norm(greatCircleU);
	normV = norm(greatCircleV);

	for (i=0; i < 3; i++) {
		greatCircleU[i] = (earthR) * greatCircleU[i] / normU;
		greatCircleV[i] = (earthR) * greatCircleV[i] / normV;
	
	}
	sunShadow = Array();
	sunShadowLatLon = Array();
	for (i=0; i < numCircle; i++) {

		theta = (i/numCircle)*Math.PI*2;
		sunShadow[i] = {x: greatCircleU[0]*Math.cos(theta) + greatCircleV[0]*Math.sin(theta), 
						y: greatCircleU[1]*Math.cos(theta) + greatCircleV[1]*Math.sin(theta),
						z: greatCircleU[2]*Math.cos(theta) + greatCircleV[2]*Math.sin(theta),
						r: 1, red: 255, green: 0, blue: 0};
		sunShadowLatLon[i] = {lon: Math.atan2(sunShadow[i].x, sunShadow[i].z),
						  lat: Math.PI/2 - Math.acos(sunShadow[i].y / earthR )};
	}

	// Velocity calculation (apprixmate)
	velocity =   (satObj[currentSat].x - pold.x) * (satObj[currentSat].x - pold.x) + (satObj[currentSat].y - pold.y) * (satObj[currentSat].y - pold.y) + (satObj[currentSat].z - pold.z) * (satObj[currentSat].z - pold.z);

	velocity = 1000 * Math.sqrt(velocity)/(tNowMs-told); // pixels/s
	velocity = metersToPix * velocity;				 // meters/s

	altitude = Math.sqrt(satObj[currentSat].x * satObj[currentSat].x + satObj[currentSat].y * satObj[currentSat].y + satObj[currentSat].z * satObj[currentSat].z);

	altitude = (altitude - earthR) * metersToPix;

	// Get visible horizon of the current satellite
	horizonPoints = visibleHorizon(satObj[currentSat], altitude);
	horizonLatLon = returnLatLon(horizonPoints);

	pold.x = satObj[currentSat].x; pold.y = satObj[currentSat].y; pold.z = satObj[currentSat].z;
	told = tNowMs;
	

	// Rotate the longitude lines to match the current time of day
	for (i=0; i < longiObj.length; i++) {
		rotateParticle(longiObj[i], 0, dayFraction * 2 * Math.PI , 0);
	}

	// Color the longitude line points according to whether they are presently lit by the sun or not.
	hemisphereDistSq = earthR * earthR + earthR * earthR;
	for (i = 0; i < longiObj.length; i ++) {
		var distSq = (xsunV - longiObj[i].x) * (xsunV - longiObj[i].x) + (ysunV - longiObj[i].y) * (ysunV- longiObj[i].y)  + (zsunV - longiObj[i].z) * (zsunV - longiObj[i].z);
		if (distSq >= hemisphereDistSq) {
			longiObj[i].red = 100; longiObj[i].green = 100; longiObj[i].blue = 100;
		} else {
			longiObj[i].red = 255; longiObj[i].green = 255; longiObj[i].blue = 255;
		}

	}

	rotateParticle(homeObj[0], 0, rotEquinox, 0); // Correctly position the user's location on the globe using the current time and the equinox info


	// Perigee Point
	//perigee.x = radius*Math.cos(3*Math.PI/2 + aoPerigee);
	//perigee.y = 0;
	//perigee.z = radius*Math.sin(3*Math.PI/2 + aoPerigee);


	//rotateParticle(perigee, 0, raAscending, inclination);

	//drawParticle(perigee.x, perigee.y, perigee.z, 3, "50, 255, 50");

	//view = view + 2*Math.PI/1080;
	xDir = new createVector([0, 0, 0],[xsunV, ysunV, zsunV,], 100, 1, [255, 0, 0]);
	
	points = Array();
	for (i=0; i < activeSats.length; i++ ) points = points.concat(satOrbit[activeSats[i]]);
	for (i=0; i < activeSats.length; i++ ) points = points.concat(satObj[activeSats[i]]);
	points = points.concat(earthObj, longiObj, homeObj, horizonPoints, xDir.pts);
	// Assign key (indices) to the objects for later re-ordering
	for (i=0; i < points.length; i++ ) {
		points[i].key = i; 
	}

	//background fill and rendering
	if (draw3D == true ) {
		//context.fillStyle = canvasBG;
		//context.fillStyle = 'rgba(0, 0, 0, 0)';
		//context.fillRect(0,0,width,height);
		context.clearRect(0, 0, width, height);
		renderObjs(context, points, view);
	}

	// Reset earth points, since we do a total rotation to the current time, not an incrimental rotation
	for (i=0; i < longiObj.length; i++) {
		rotateParticle(longiObj[i], 0, -dayFraction * 2 * Math.PI , 0);
	}

	rotateParticle(homeObj[0], 0, -rotEquinox, 0);

	// Get numerical satellite information
	// 	Right Ascension (RA) is always measured Eastward (negative direction in our coordinate system)
	// the RA and DEC calculations need to be updated for the more general elliptical case
	// var satRA = (raAscending - Math.PI/2 + aoPerigee + totalAnomoly) / (2*Math.PI) - Math.floor( (raAscending - Math.PI/2 + aoPerigee + totalAnomoly) / (2*Math.PI) ); // This is an easy calculation because we aligned the vernal equinox to our x-axis
	// var satDEC = Math.atan( Math.sqrt( p[sat].x*p[sat].x + p[sat].z*p[sat].z )/p[sat].y) > 0 ? Math.PI/2 - Math.atan( Math.sqrt( p[sat].x*p[sat].x + p[sat].z*p[sat].z )/p[sat].y) : -Math.PI/2 - Math.atan( Math.sqrt( p[sat].x*p[sat].x + p[sat].z*p[sat].z )/p[sat].y);


	//elementName.innerHTML 		= objectName;
	elementEpoch.innerHTML  	= epochObj.toUTCString();
	elementTime.innerHTML   	= tnow.toUTCString();
	elementEccentric.innerHTML 	= (Math.round(E*1000000)/1000000).toFixed(6) + " rads";

	if (currentSat == defaultSat) {
		elementISSspeed.innerHTML	= ((velocity/1000) * 60 *60 * 0.621371).toFixed(2) + ' mph'; // Conversion from m/s to mph
		elementISSalt.innerHTML		= ((altitude/1000) * 0.621371).toFixed(2) + ' miles';			// Conversion from m to miles
	} else {
		elementISSspeed.innerHTML	= '17,000 mph';
		elementISSalt.innerHTML		= '256 miles';
	}

	elementVelocity.innerHTML	= ((velocity/1000)).toFixed(4) + " km/s";
	elementAltitude.innerHTML 	= (altitude/1000).toFixed(4) + " km";
	//elementSatRA.innerHTML  = Math.round(360 * satRA *1000000)/1000000;
	//elementSatDEC.innerHTML = Math.round( (180/(Math.PI)) * satDEC *1000000)/1000000;
	//console.log("Right Ascension: "+ satRA * 360)
	//console.log( Math.floor((24*(satRA/(2*Math.PI)))) + "h, " + Math.floor( 60*((24*(satRA/(2*Math.PI)))%24) ) + "m, " + Math.floor( 60*((24*60*(satRA/(2*Math.PI)))%(24 * 60)) ) + "s");

}

function readTLE(tleStr) {
	var lines = tleStr.split('\n');
	var TLE1 = lines[0].split(/\s+/);
	var TLE2 = lines[1].split(/\s+/);
	var TLE3 = lines[2].split(/\s+/);

	this.name 		 = lines[0];
	this.year 		 = parseInt('20'+TLE2[3].slice(0,2));
	this.epoch 		 = parseFloat(TLE2[3].slice(2));
	this.inclination = parseFloat(TLE3[2]);
	this.raAscending = parseFloat(TLE3[3]);
	this.eccentricity = parseFloat('0.'+TLE3[4]);
	this.aoPerigee	 = parseFloat(TLE3[5]);
	this.meanAnomoly = parseFloat(TLE3[6]);
	this.meanMotion  = parseFloat(TLE3[7]);

}

function constructTLEObjs(dataStr) {
	var objArray = Array();
	var lines = dataStr.split('\n'); 

	imax = Math.floor(lines.length/3);
	for (i = 0; i < imax; i++) {
		tleStr = lines[i*3] + '\n' + lines[i*3+1] + '\n' + lines[i*3+2];
		objArray[i] = new readTLE(tleStr);
	}

	return objArray;
}


//function update() {
//	satellite();
//	if ( draw3D == true) { timer = setTimeout(update, tUpdate) };
//}

var canvasHolder = document.getElementById("canvasHolder");

if ( draw3D == true ) {

	var canvasObj = document.createElement('canvas');
	canvasObj.id = "satelliteCanvas";
	canvasObj.width = 200;
	canvasObj.height = 200;

	canvasHolder.appendChild(canvasObj);

	//canvasObj canvasObj = document.getElementById("satelliteCanvas");
	var context = canvasObj.getContext("2d");


	var height = canvasObj.height;
	var width =  canvasObj.width;
} else {
	height = 400;
	width = 400;
}

var scale = height/900; 

var nSats = 0;
var satOrbit = Array();
var satObj = Array();
var perigee = new Object();
var rotEquinox;
var dragging = false;
var mx;
var my;
var E;
var colorStyle = 'dark';
var canvasBG;

var G 			= 6.67384 * Math.pow(10,-11);	// m^3 kg^-1 s^-2
var earthMass 	= 5.97219 * Math.pow(10, 24);	// kg
var earthRadius = 6371;							// km

var numCircle = Math.ceil(scale*800);
var refRadius = scale*300; 		// pixels, everything gets scaled relative to the semi-major axis of the orbit

var numLongi = Math.ceil(scale*200);

var xc = height/2;
var yc = width/2;

var longiCircle = Array(numLongi);
var longiObj   = Array();

// Get user's coordinates
var homeLat		= elementULat.value;
var homeLongi	= elementULongi.value;

getUsersCoordinates();

var homeObj		= Array();
var satInfo		= Array();
var activeSats  = Array();
var maxDepth    = 0;


// Define some global vars for a finite-difference velocity calculation
var told;
var tnew;
var pold = new Object;
var pnew = new Object;
var altitude;
var velocity;

// Define a reference state to orient earth correctly in the equitorial coordinate system
//	One way to do this is to use the vernal equinox. At this day and time, the sun passes
//	through the ascending node where the elliptical and celestial equator intersect.
//	
//	Knowing the date and time of this event, and that UTC is measured at 0 degrees longitude,
//	it is a straightforward, albeit imperfect, calculation to find the longitude of the
//	vernal equinox. We can make all further calculations relative to this point
//
//	Right Ascension begins at the vernal equinox 

// Wednesday March 20th, 11:02 UTC = Vernal Equinox
// March 20th, 2014, 16:57 UTC = Vernal Equinox
var vernalEquinox = Date.UTC(2014, 2, 20, 16, 57, 0, 0, 0);  			  // UTC time (0 degrees longitude) when Vernal Equinox occurs
var equinoxNoon = Date.UTC(2014, 2, 20, 12, 0, 0, 0, 0)                   // UTC noon (when sun is directly overhead)
var equinoxLongi  = 2 * Math.PI * (equinoxNoon - vernalEquinox)/dayMs; // Difference between noon, and time of equinox, divided by time in day, multiplied by radians of one rotation
                                                                          //   gives the longitude where the vector from the center of Earth points to the sun at the time of equinox
var TLEelement   = document.getElementById("TLEData");
var dropdownMenu = document.getElementById("selectSatellites");

var TLEdata = TLEelement.innerHTML;
var TLEobjs = constructTLEObjs(TLEdata);

for (i = 0; i < TLEobjs.length; i++) {
	var opt = TLEobjs[i];
	var el = document.createElement("option");
	if (i == 0) {
		el.selected = "selected";
	}
	el.textContent = opt.name;
	el.value 	   = i;
	dropdownMenu.appendChild(el);

}
var defaultSat = 0;
var currentSat = defaultSat;
var defaultObjectName = TLEobjs[defaultSat].name;

var objectName = defaultObjectName;


satInfo[defaultSat] = new setSatellite(TLEobjs[defaultSat]);
activeSats.push(defaultSat);

var metersToPix = satInfo[defaultSat].a /refRadius; 	//meters/pixel 

var view = {x: 0*2 * Math.PI * homeLat/360, y: 0*2 * Math.PI * homeLongi/360 - equinoxLongi }

init();



if ( draw3D == true ) {
	//fixes a problem where double clicking causes text to get selected on the canvas
	canvasObj.addEventListener('selectstart', function(e) { e.preventDefault(); return false; }, false);

	//Prevent default touch events on the canvas
	canvasObj.addEventListener("touchmove", function(e) { e.preventDefault(); return false; }, false);

	canvasObj.addEventListener('mousedown', function(e) {
		mx = e.pageX;
		my = e.pageY; 
		dragging = true; 
		tUpdate = 30;	// for smooth drawing while moving the object
		clearTimeout(timer)
		timer = setTimeout(update, tUpdate);
	}, true);

	canvasObj.addEventListener('mousemove', function(e) {
		
		if (dragging) {
			view.y = view.y + ((e.pageX - mx) / (2 *width) ) * Math.PI * 2;
			view.x = view.x + ((e.pageY - my) / (2 *width) ) * Math.PI * 2;
			//console.log(e.pageX)
			mx = e.pageX;
			my = e.pageY;

		}

	}, true);

	canvasObj.addEventListener('mouseup', function(e) {
		dragging = false; 
		tUpdate = 1000;
		clearTimeout(timer)
		timer = setTimeout(update, tUpdate);
	}, true)

	canvasObj.addEventListener('dblclick', function(e) {
		view.x = view.y = 0;
		dragging = false;
		tUpdate = 1000;
		clearTimeout(timer)
		timer = setTimeout(update, tUpdate);
	}, true)
}

//if (draw3D) { var timer = setTimeout(update, tUpdate) };

//Draw satellite label

	//padding = 4*scale;
	//spacing = 20*scale;
	//idFontSize = 12;

	//context.textAlign = "left";
	//context.textBaseline = "middle";
	//context.font = idFontSize + "px sans-serif";
	//context.fillText("ISS", xc + tempObjs[sat].x + spacing, yc -p[sat].y);

	//context.strokeStyle = "rgba(255, 255, 255, 1)";
	//context.beginPath();
	//context.moveTo(xc + p[sat].x + spacing + p[sat].r - padding*2 , yc -p[sat].y);
	//context.lineTo(xc + p[sat].x + p[sat].r + padding, yc -p[sat].y );
	//context.stroke();


