#version 330 core
precision highp float;
in vec2 coord;
out vec4 fragColor;

uniform vec3 cameraPosition;
uniform mat4 projectiveInverse;
uniform float fTime;
uniform float magicNumber;
uniform float magicNumberLin;
uniform float mandelboxScale;
uniform float kaleido;
uniform int deType;
uniform vec3 bassHole;
uniform vec3 midsHole;
uniform vec3 highHole;
uniform vec3 avgMidsHole;
uniform vec3 avgHighHole;
uniform vec2 avgMidsPlot;
uniform vec2 avgHighPlot;

const float pi = 3.14159265358;
const float tau = 2.0*pi;
const float e = 2.718281828;
const int maxIterations = 150;
const float epsilon = 0.00001;
const vec3 dirX = vec3(1.0, 0.0, 0.0);
const vec3 dirY = vec3(0.0, 1.0, 0.0);
const vec3 dirZ = vec3(0.0, 0.0, 1.0);

// Phong lighting
const float ambientStrength = 0.35;
const vec3 lightDir = normalize(vec3(1.0, -0.6, -1.5));
const vec3 lightColor = vec3(0.95, 0.95, 0.95);
const vec3 ambientLight = ambientStrength * lightColor;

mat3 buildRot3(vec3 u, float theta)
{
	float c = cos(theta);
	float cC = 1.0 - c;
	float s = sin(theta);
	float sC = 1.0 - s;
	return mat3(
		c+u.x*u.x*cC, u.y*u.x*cC+u.z*s, u.z*u.x*cC-u.y*s,
		u.x*u.y*cC-u.z*s, c+u.y*u.y*cC, u.z*u.y*cC+u.x*s,
		u.x*u.z*cC+u.y*s, u.y*u.z*cC-u.x*s, c+u.z*u.z*cC
	);
}

vec2 powc(float x, vec2 s)
{
	// x^s = e^(log(x)s) = e^(log(x)Re(s) + log(x)Im(s)i) = x^Re(s) * (cos(log(x)Im(s)) + sin(log(x)Im(s))i)
	float theta = log(x) * s.y;
	return pow(x, s.x) * vec2(cos(theta), sin(theta));
}

vec2 Normalize(vec2 v) {
	float r2 = v.x*v.x + v.y*v.y;
	if (r2 < 0.000001) {
		return v;
	}
	return normalize(v);
}

float getAngle(vec2 s)
{
	float theta = 0.0;
	if (s.y < 0.0) {
		s *= -1.0;
		theta = pi;
	}

	s = Normalize(s);
	if (s.x >= 0.0)
		return theta + asin(s.y);
	else
		return theta + pi - asin(s.y);
}

vec2 multInv(vec2 s)
{
	return vec2(s.x, -s.y)/dot(s, s);
}

vec2 multC(vec2 a, vec2 b)
{
	return vec2(a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x);
}

vec2 sinc(vec2 z)
{
	// z = a + bi
	float sina = sin(z.x);
	float cosa = cos(z.x);
	return 0.5*(exp(z.y)*vec2(sina, cosa) - exp(-z.y)*vec2(-sina, cosa));
}

vec2 expc(vec2 z){return powc(e, z);}

vec2 cpowr(vec2 c, float x)
{
	float r2 = c.x*c.x + c.y*c.y;
	float theta = x * getAngle(c);
	return pow(r2, x/2.0) * vec2(cos(theta), sin(theta));
}

float bound(float x, float b) {
	return mod(x + b, 2.0*b) - b;
}
float boundReflect(float x, float b) {
	float r = mod(x + b, 4.0*b);
	if (r < 2.0*b) {
		return r - b;
	} else {
		return 3.0*b - r;
	}
}

vec3 gradient;
vec4 orbitTrap;
const vec3 cutNorm = normalize(vec3(-1.0, 1.0, -1.0));
float distanceEstimator(vec3 t)
{
	orbitTrap = vec4(1.0, 1.0, 1.0, 1.0);

	//*/
	// Mandelbox
	if (deType == 1) {
		const int maxIterations = 6;
		vec3 s = t;
		float DEfactor = 1.0;
		float r2 = 1.0;
		const float maxR2 = 12.0;
		const float BVR = sqrt(maxR2);
		for (int i = 0; i < maxIterations; i++) {
			if(s.x>1.0){s.x=2.0-s.x;}else if(s.x<-1.0){s.x=-2.0-s.x;}
			if(s.y>1.0){s.y=2.0-s.y;}else if(s.y<-1.0){s.y=-2.0-s.y;}
			if(s.z>1.0){s.z=2.0-s.z;}else if(s.z<-1.0){s.z=-2.0-s.z;}

			r2 = dot(s, s);
			if (r2 < 0.25) {
				s *= 4.0;
				DEfactor *= 4.0;
			} else if(r2 < 1.0) {
				s /= r2;
				DEfactor /= r2;
			}

			orbitTrap.x = min(orbitTrap.x, length(s/BVR - bassHole)/2.0);
			orbitTrap.y = min(orbitTrap.y, length(s/BVR - midsHole)/2.0);
			orbitTrap.z = min(orbitTrap.z, length(s/BVR - highHole)/2.0);
			//orbitTrap.w = min(orbitTrap.w, length(s));

			s = s*mandelboxScale + t;
			DEfactor = DEfactor*abs(mandelboxScale) + 1.0;
		
			if(r2 > maxR2) break;
		}
		return (length(s)-BVR)/abs(DEfactor);
	}//*/
	//*/
	// Mandelbulb
	else if (deType == 2) {
		const int maxIterations = 4;
		const float reScale = 0.85;
		t = vec3(boundReflect(t.x, 9.0), boundReflect(t.y, 9.0), boundReflect(t.z, 9.0));
		t *= reScale;
		vec3 s = t;
		float power = 6.0 + 9.0*magicNumberLin;
		float dr = 1.0;
		float r = 0.0;
		for (int i = 0; i < maxIterations; i++) {
			r = length(s);
			const float b = 1.5;
			if (r > b) break;

			float theta = acos(s.z/r);
			float phi = atan(s.y, s.x);
			dr = pow(r, power-1.0)*power*dr + 1.0;

			r = pow(r, power);
			theta *= power;
			phi *= power;

			s = r*vec3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
			s += t;

			orbitTrap.xyz = min(orbitTrap.xyz, abs(s - (highHole + midsHole)/1.75));
		}
		return min(0.5*log(r)*r/dr / reScale, 4.5);
	}//*/
    //*/
    // Knighty's Pseudo Klienian*
	else if (deType == 3) {
		const int maxIterations = 4;
		const float reScale = 0.5;
		vec3 s = reScale*t;

		float anim = 1.25 + 0.08*magicNumber;
		float scale = 1.0;
		float theta = 0.25*fTime;
		mat2 rotato = mat2(cos(theta), sin(theta), -sin(theta), cos(theta));

		for(int i = 0; i < maxIterations; i++)
		{
			if (i == 2) {
				s.xz *= rotato;
			}

			s = -1.0 + 2.0*fract(0.5*s + 0.5);

			float r2 = dot(s,s);
		
			float k = anim/r2;
			s *= k;
			scale *= k;

			orbitTrap.x = min(orbitTrap.x, length(s/1.5 - bassHole)/2.0);
			orbitTrap.y = min(orbitTrap.y, length(s/1.5 - midsHole)/2.0);
			orbitTrap.z = min(orbitTrap.z, length(s/1.5 - highHole)/2.0);
		}
	
		return (0.25*abs(s.y)/scale) / reScale;
	} //*/
	//*/
	// Menger
	else if (deType == 4) {
		const int maxIterations = 5;

		const float reScale = 0.425;
		vec3 s = t*reScale;

		s = vec3(boundReflect(s.x, 4.0), boundReflect(s.y, 4.0), boundReflect(s.z, 4.0));

		s = s + 0.5; //center it by changing position and scale
		float xx=abs(s.x-0.5)-0.5, yy=abs(s.y-0.5)-0.5, zz=abs(s.z-0.5)-0.5;
		float d1=max(xx,max(yy,zz)); //distance to the box
		float d=d1; //current computed distance
		float p=1.0;
		//float mengerScale = 2.55 + 0.85*magicNumberLin;
		float mengerScale = 3.0;
		float halfScale = mengerScale / 2.0;

		orbitTrap.xyz = abs(vec3(xx/1.25, yy/1.25, zz/1.25));

		float theta = 0.04*fTime;
		//mat3 rotato = buildRot3(normalize(bassHole), theta);
		mat3 rotato = buildRot3(normalize(vec3(1, 2, 2)), theta);

		for (int i = 0; i < maxIterations; i++) {
			p *= mengerScale;
			float xa = mod(s.x*p, mengerScale);
			float ya = mod(s.y*p, mengerScale);
			float za = mod(s.z*p, mengerScale);

			float xx=0.5-abs(xa-halfScale), yy=0.5-abs(ya-halfScale), zz=0.5-abs(za-halfScale);
			d1=min(max(xx,zz),min(max(xx,yy),max(yy,zz))) / p; //distance inside the 3 axis-aligned square tubes

			d=max(d,d1); //intersection

			if (i % 2 == 0) {
				const float rat = 0.68;
				vec3 q = vec3(xx, yy, zz)/mengerScale;
				vec3 col = vec3(length(q - bassHole)/3.75, length(q - midsHole)/3.75, length(q - highHole)/3.75);
				orbitTrap.xyz = rat*orbitTrap.xyz + (1.0 - rat)*col;
			}

			s *= rotato;
		}
		return d/reScale;
	}//*/
	//*/
	// Knighty's Kaleidoscopic IFS
	else if (deType == 5) {
		const int maxIterations = 8;
		const float reScale = 0.65;
		const float scale = 2.0;
		t = vec3(boundReflect(t.x, 9.0), boundReflect(t.y, 9.0), boundReflect(t.z, 9.0));
		vec3 s = reScale*t;
		const vec3 center = vec3(sqrt(0.5), sqrt(0.3), sqrt(0.2));
		float r2 = dot(s, s);
		float DEfactor = 1.0;

		float theta = 0.5*fTime;
		mat3 rotato1 = buildRot3(normalize(avgMidsHole), theta);
		theta = 0.1105*pi*boundReflect(6.0*magicNumberLin, 1.0);
		mat3 rotato2 = buildRot3(normalize(avgHighHole), theta);

		for(int i = 0; i < maxIterations && r2 < 1000.0; i++) {
			s *= rotato1;

			if(s.x+s.y<0.0){float x1=-s.y;s.y=-s.x;s.x=x1;}
			if(s.x+s.z<0.0){float x1=-s.z;s.z=-s.x;s.x=x1;}
			if(s.y+s.z<0.0){float y1=-s.z;s.z=-s.y;s.y=y1;}

			s *= rotato2;

			s = scale*s - (scale - 1.0)*center;
			r2 = dot(s, s);

			orbitTrap.x = min(orbitTrap.x, length(s - bassHole)/2.0);
			orbitTrap.y = min(orbitTrap.y, length(s - midsHole)/2.0);
			orbitTrap.z = min(orbitTrap.z, length(s - highHole)/2.0);

			DEfactor *= scale;
		}
		return (sqrt(r2) - 2.0) / DEfactor / reScale;
	}//*/
	else {
		return 1000.0;
	}
}

const float maxBrightness = 1.585;
const float maxBrightnessR2 = maxBrightness*maxBrightness;
vec4 scaleColor(float si, vec3 col) {
	col *= pow(1.0 - si/float(maxIterations), 1.1);
	if(dot(col, col) > maxBrightnessR2) {
		col = maxBrightness*normalize(col);
	}
	return vec4(col, 1.0);
}

vec4 phongLighting(vec3 c) {
	vec3 diffuse = max(dot(normalize(gradient), -lightDir), 0.0) * lightColor;
	return vec4((ambientLight + diffuse) * c, 1.0);
}

void main(void)
{
	const float d = 1.0;
	const float near = 1.0;
	const float far = 2.0;
	const float projectionConstant = d*(far+near)/(far-near) - (2.0*far*near)/(far-near);
	const float maxDistance = 100.0;
	//vec2 tempCoord = coord-avgHighPlot-avgMidsPlot;
	vec2 tempCoord = coord;
	float magicTheta = boundReflect(getAngle(tempCoord),(1.0-kaleido)*tau + kaleido*pi/6.0);
	vec2 tCoord = length(tempCoord) * vec2(cos(magicTheta), sin(magicTheta));
	magicTheta = pi*magicNumber;
	float choiceDot = dot(coord.xy, vec2(cos(magicTheta), sin(magicTheta)));
	vec2 choice = choiceDot > 0.0 ? avgHighPlot : avgMidsPlot;
	float otherThing = 2.2*dot(choice, tCoord.xy);
	vec2 newCoord = tCoord.xy + pow(abs(choiceDot), 1.145)*otherThing*choice;
	newCoord = vec2(boundReflect(newCoord.x, 1.0), boundReflect(newCoord.y, 1.0));
	vec3 direction = normalize((projectiveInverse * vec4(newCoord.x*d, newCoord.y*d, projectionConstant, d)).xyz);
	const float minTravel = 0.65;
	vec3 pos = deType == 2 ? cameraPosition : cameraPosition + minTravel*direction;
	const float hitDistance = 0.00001;
	float minDist = maxDistance;
	float minI = 0.0;
	float travel = minTravel;
	magicTheta = 2.0*getAngle(newCoord);
	vec2 newestCoord = length(newCoord)*vec2(cos(magicTheta), sin(magicTheta));
	float newMagic = (magicNumber + 2.0)/2.0 + 1.75*dot(coord, avgMidsPlot);
	for (int i = 0; i < maxIterations; i++) {
		float dist = newMagic - travel;
		if(dist > 0.0) {
			pos += dist*direction;
			travel += dist;
			i++;
		}
		dist = distanceEstimator(pos);

		if (dist <= hitDistance) {
			float smoothI = float(i);// - (log(travel)/10.0);
			orbitTrap.x = pow(orbitTrap.x, 0.725);
			orbitTrap.y = pow(orbitTrap.y, 0.725);
			orbitTrap.z = pow(orbitTrap.z, 0.725);
			fragColor = scaleColor(smoothI, orbitTrap.xyz);
			//fragColor = phongLighting(scaleColor(smoothI, orbitTrap.xyz).xyz);
			return;
		}

		if (dist < minDist) {
			minDist = dist;
			minI = float(i);
		}

		pos += 0.99*dist*direction;
		travel += dist;
		if (travel > maxDistance) {
			direction = normalize((projectiveInverse * vec4(coord.x*d, coord.y*d, projectionConstant, d)).xyz);
			vec3 sinDir = sin(100.0*direction);
			vec3 base = vec3(exp(-3.0*length(sin(pi * bassHole + 1.0) - sinDir)), exp(-4.0*length(sin(e * midsHole + 1.3) - sinDir)), exp(-3.0*length(sin(9.6*highHole + 117.69420) - sinDir)));
			fragColor = vec4((deType == 0 ? 0.75 : 0.45) * base, 1.0);
			return;
		}
	}
	fragColor = vec4(0.0, 0.0, 0.0, 1.0);
	return;
}