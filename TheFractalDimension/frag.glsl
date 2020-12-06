#version 330 core
precision highp float;
in vec2 coord;
out vec4 fragColor;

uniform vec3 cameraPosition;
uniform mat4 projectiveInverse;
uniform float magicNumber;
uniform float mandelboxScale;
uniform bool question;
uniform vec3 bassHole;
uniform vec3 midsHole;
uniform vec3 highHole;
uniform vec2 avgMidsHole;
uniform vec2 avgHighHole;

const float pi = 3.14159265358;
const float e = 2.718281828;
const int maxIterations = 40;
const float epsilon = 0.00001;
const vec3 dirX = vec3(1.0, 0.0, 0.0);
const vec3 dirY = vec3(0.0, 1.0, 0.0);
const vec3 dirZ = vec3(0.0, 0.0, 1.0);

// Phong lighting
const float ambientStrength = 0.35;
const vec3 lightDir = normalize(vec3(1.0, -0.6, -1.5));
const vec3 lightColor = vec3(0.95, 0.95, 0.95);
const vec3 ambientLight = ambientStrength * lightColor;

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

const float power = 9.0;
vec3 gradient;
vec4 orbitTrap;
const vec3 cutNorm = normalize(vec3(-1.0, 1.0, -1.0));
float potential(vec3 t)
{
	vec3 s = t;
	const int maxIterations = 20;
	for(int i = 1; i < maxIterations; i++)
	{
		return log(length(s))/pow(power, float(i));
	}
	return 0.0;
}
float d_basic(vec3 t)
{
	const float thick = 0.0015;
    float r = length(t.x);
    return max(r, -(length(t.xy)*t.z - thick)/sqrt(dot(t, t) + thick));
}
float distanceEstimator(vec3 t)
{
	//t += 0.25*t*vec3(exp(cos(dot(bassHole,bassHole)) + 0.5), exp(cos(dot(midsHole,midsHole)) + 0.5), exp(cos(dot(highHole,highHole)) + 0.5));
	//t = vec3(bound(t.x, 4.0), bound(t.y, 4.0), bound(t.z, 4.0));
	vec3 s = t;
    //vec3 s = vec3(bound(t.x, 4.0), bound(t.y, 4.0), bound(t.z, 4.0));
	
	const int maxIterations = 7;
	orbitTrap = vec4(1.0, 1.0, 1.0, 1.0);

	// Balls
	//return length(t) - 0.1;

	/*/
	float pot = potential(t);
	if(pot <= 0.0) return 0.0;
	gradient = (vec3(potential(t + epsilon*dirX), potential(t + epsilon*dirY), potential(t + epsilon*dirZ)) - pot)/epsilon;
	return (0.5/exp(pot))*sinh(pot)/length(gradient);//*/

	/*/
	// Mandelbulb OLD
	float dr = 1.0;
	float r = 0.0;
	for (int i = 0; i < maxIterations; i++) {
		r = length(s);
		if (r > 1.25) break;

		float theta = acos(s.z/r);
		float phi = atan(s.y, s.x);
		dr = pow(r, power-1.0)*power*dr + 1.0;

		float zr = pow(r, power);
		theta *= power;
		phi *= power;

		s = zr*vec3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
		s += t;

		orbitTrap.x = min(orbitTrap.x, length(s - dirX));
		orbitTrap.y = min(orbitTrap.y, length(s - dirY));
		orbitTrap.z = min(orbitTrap.z, length(s - dirZ));
	}
	return 0.5*log(r)*r/dr;//*/

    /*/
    // Knighty's Pseudo Klienian*
    const vec3 cellSize = vec3(0.808001, 0.808, 1.167);
	float DEfactor = 1.0;
	for(int i = 0; i < maxIterations; i++) {
		s = 2.0 * clamp(s, -cellSize, cellSize) - s;
		//Inversion
		float r2 = dot(s, s);
		float k = max(1.0 / r2, 1.0);
		s *= k; DEfactor *= k;

		orbitTrap.x = min(orbitTrap.x, length(s - dirX));
		orbitTrap.y = min(orbitTrap.y, length(s - dirY));
		orbitTrap.z = min(orbitTrap.z, length(s - dirZ));
	}
	return (0.363001*d_basic(s)/abs(DEfactor));//*/

	/*/
	// Menger
	s = s/2.0 + 0.5; //center it by changing position and scale
	float xx=abs(s.x-0.5)-0.5, yy=abs(s.y-0.5)-0.5, zz=abs(s.z-0.5)-0.5;
	float d1=max(xx,max(yy,zz)); //distance to the box
	float d=d1; //current computed distance
	float p=1.0;
	for (int i = 0; i < maxIterations; i++) {
		float xa = mod(3.0*s.x*p,3.0);
		float ya = mod(3.0*s.y*p,3.0);
		float za = mod(3.0*s.z*p,3.0);
		p*=3.0;

		float xx=0.5-abs(xa-1.5), yy=0.5-abs(ya-1.5), zz=0.5-abs(za-1.5);
		d1=min(max(xx,zz),min(max(xx,yy),max(yy,zz))) / p; //distance inside the 3 axis-aligned square tubes

		d=max(d,d1); //intersection
	}
	return max(d, 0.0*dot(cutNorm, t));//*/

    //*/
    // Mandelbox
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

		orbitTrap.x = min(orbitTrap.x, length(s - bassHole)/(2.0 * BVR));
		orbitTrap.y = min(orbitTrap.y, length(s - midsHole)/(2.0 * BVR));
		orbitTrap.z = min(orbitTrap.z, length(s - highHole)/(2.0 * BVR));
		//orbitTrap.x = min(orbitTrap.x, length(s + dirX));
		//orbitTrap.y = min(orbitTrap.y, length(s + dirY));
		//orbitTrap.z = min(orbitTrap.z, length(s + dirZ));
		//orbitTrap.x = min(orbitTrap.x, length(s - BVR*dirX)/BVR);
		//orbitTrap.y = min(orbitTrap.y, length(s - BVR*dirY)/BVR);
		//orbitTrap.z = min(orbitTrap.z, length(s - BVR*dirZ)/BVR);
		//orbitTrap.w = min(orbitTrap.w, length(s));

		s = s*mandelboxScale + t;
		DEfactor = DEfactor*abs(mandelboxScale) + 1.0;
		
		if(r2 > maxR2) break;
	}
	return (length(s)-BVR)/abs(DEfactor); //*/
}

const float maxBrightness = 1.4;
const float maxBrightnessR2 = maxBrightness*maxBrightness;
vec4 scaleColor(float si, vec3 col) {
	col *= pow(1.0 - si/float(maxIterations), 0.275);
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
	const float maxDistance = 4.0;
	float magicTheta = pi*magicNumber;
	float choiceDot = dot(coord.xy, vec2(cos(magicTheta), sin(magicTheta)));
	vec2 choice = choiceDot > 0.0 ? avgHighHole : avgMidsHole;
	float otherThing = 2.2*dot(choice, coord.xy);
	vec2 newCoord = coord.xy + pow(abs(choiceDot), 1.145)*otherThing*choice;
	newCoord = vec2(boundReflect(newCoord.x, 1.0), boundReflect(newCoord.y, 1.0));
	vec3 direction = normalize((projectiveInverse * vec4(newCoord.x*d, newCoord.y*d, projectionConstant, d)).xyz);
	const float minTravel = 0.63;
	vec3 pos = cameraPosition + minTravel*direction;
	const float hitDistance = 0.00001;
	float minDist = maxDistance;
	float minI = 0.0;
	float travel = minTravel;
	magicTheta = 2.0*getAngle(newCoord);
	vec2 newestCoord = length(newCoord)*vec2(cos(magicTheta), sin(magicTheta));
	float newMagic = (magicNumber + 2.0)/2.0 + 1.75*dot(newestCoord, avgMidsHole);
	for (int i = 0; i < maxIterations; i++) {
		float dist = question ? 1000.0 : newMagic - travel;
		if(dist > 0.0) {
			pos += dist*direction;
			travel += dist;
			i++;
		}
		dist = min(0.99 * distanceEstimator(pos), 1.0);

		if (dist <= hitDistance) {
            float smoothI = float(i);// - (log(travel)/10.0);
			orbitTrap.x = pow(orbitTrap.x, 0.725);
			orbitTrap.y = pow(orbitTrap.y, 0.725);
			orbitTrap.z = pow(orbitTrap.z, 0.725);
			//orbitTrap = vec4(mix(orbitTrap.xyz, abs(pos), exp(-(pos.x*pos.x + pos.y*pos.y + pos.z*pos.z))), orbitTrap.w);
			//orbitTrap = vec4(mix(orbitTrap.xyz, abs(pos), exp(-(abs(pos.x) + abs(pos.y) + abs(pos.z)))), orbitTrap.w);
			fragColor = scaleColor(smoothI, orbitTrap.xyz);
			return;
		}

		if (dist < minDist) {
			minDist = dist;
			minI = float(i);
		}

		pos += dist*direction;
		travel += dist;
		if (travel > maxDistance) {
			direction = normalize((projectiveInverse * vec4(coord.x*d, coord.y*d, projectionConstant, d)).xyz);
			vec3 sinDir = sin(100.0*direction);
			vec3 base = vec3(exp(-3.0*length(sin(pi * bassHole + 1.0) - sinDir)), exp(-4.0*length(sin(e * midsHole + 1.3) - sinDir)), exp(-3.0*length(sin(9.6*highHole + 117.69420) - sinDir)));
			fragColor = vec4((question ? 0.65 : 0.4) * base, 1.0);
			return;
		}
	}
	fragColor = vec4(0.0, 0.0, 0.0, 1.0);
	return;
}