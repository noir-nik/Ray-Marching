#version 330

uniform vec2 iResolution;
out vec4 fragColor;

uniform int uMaxSteps;
uniform float uMaxDist;
uniform float uEps;
uniform int uAa;
uniform int uReflect;

// Camera
uniform vec3 uCamPos;
uniform vec3 uLookAt;
uniform float uFov;

// Ops
vec2 opUnion(vec2 d1, vec2 d2) {
	return (d1.x < d2.x) ? d1 : d2;
}
vec2 opSubtraction(vec2 d1, vec2 d2) {
	return max(-d1, d2);
}
vec2 opIntersection(vec2 d1, vec2 d2) {
	return (d1.x > d2.x) ? d1 : d2;
}
vec2 opXor(vec2 d1, vec2 d2) {
	return max(min(d1, d2), -max(d1, d2));
}

// SMOOTH BOOLEAN OPERATORS //

// Smooth Union
// d1: signed distance to shape 1
// d2: signed distance to shape 2
// k: smoothness value for the trasition
float opUS(float d1f, float d2f, float kf) {
	float h = clamp(0.5 + 0.5 * (d2f - d1f) / kf, 0.0, 1.0);
	float dist = mix(d2f, d1f, h) - kf * h * (1.0 - h);

	return dist;
}

// Smooth Subtraction
// d1: signed distance to shape 1
// d2: signed distance to shape 2
// k: smoothness value for the trasition
float opSS(float d1, float d2, float k) {
	float h = clamp(0.5 - 0.5 * (d2 + d1) / k, 0.0, 1.0);
	return mix(d2, -d1, h) + k * h * (1.0 - h);
}

// Smooth Intersection
// d1: signed distance to shape 1
// d2: signed distance to shape 2
// k: smoothness value for the trasition
vec4 opIS(vec4 d1, vec4 d2, float kf) {
	float h = clamp(0.5 - 0.5 * (d2.w - d1.w) / kf, 0.0, 1.0);
	float dist = mix(d2.w, d1.w, h) + kf * h * (1.0 - h);
	vec3 color = mix(d2.rgb, d1.rgb, h);

	return vec4(color.rgb, dist);
}

// SDFs
float sdSphere(vec3 p, float s) {
	return length(p) - s;
}

// generic ellipsoid - improved approximated distance
float sdEllipsoid( in vec3 p, in vec3 r ) 
{
    float k0 = length(p/r);
    float k1 = length(p/(r*r));
    return k0*(k0-1.0)/k1;
}

float sdVerticalCapsule( vec3 p, float h, float r )
{
  p.y -= clamp( p.y, 0.0, h );
  return length( p ) - r;
}

float sdBox(vec3 p, vec3 b) {
	vec3 q = abs(p) - b;
	return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);
}

// Round Box - exact
float sdRoundBox(vec3 p, vec3 b, float r) {
	vec3 q = abs(p) - b + r;
	return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0) - r;
}

// Plane - exact
float sdPlane(vec3 p, vec3 n, float h) {
  // n must be normalized
	return dot(p, n) + h;
}

// Vertical Capped Cylinder - exact
float sdCappedCylinder(vec3 p, float h, float r) {
	vec2 d = abs(vec2(length(p.xz), p.y)) - vec2(r, h);
	return min(max(d.x, d.y), 0.0) + length(max(d, 0.0));
}

// Rounded Cylinder - exact
float sdRoundedCylinder(vec3 p, float ra, float rb, float h) {
	vec2 d = vec2(length(p.xz) - 2.0 * ra + rb, abs(p.y) - h);
	return min(max(d.x, d.y), 0.0) + length(max(d, 0.0)) - rb;
}

// Cone - bound (not exact!)
float sdCone(vec3 p, vec2 c, float h) {
	float q = length(p.xz);
	return max(dot(c.xy, vec2(q, p.y)), -h - p.y);
}

// Torus - exact
float sdTorus(vec3 p, vec2 t) {
	vec2 q = vec2(length(p.xz) - t.x, p.y);
	return length(q) - t.y;
}

// Capped Torus
float sdCappedTorus( vec3 p, vec2 sc, float ra, float rb)
{
  p.x = abs(p.x);
  float k = (sc.y*p.x>sc.x*p.y) ? dot(p.xy,sc) : length(p.xy);
  return sqrt( dot(p,p) + ra*ra - 2.0*ra*k ) - rb;
}

// Gem
float sdSolidAngle(vec3 p, vec2 c, float ra) {
  // c is the sin/cos of the angle
	vec2 q = vec2(length(p.xz), p.y);
	float l = length(q) - ra;
	float m = length(q - c * clamp(dot(q, c), 0.0, ra));
	return max(l, m * sign(c.y * q.x - c.x * q.y));
}

float sdCross(vec3 p) {
	float da = sdBox(p.xyz, vec3(10.0, 1.0, 1.0));
	float db = sdBox(p.yzx, vec3(1.0, 10.0, 1.0));
	float dc = sdBox(p.zxy, vec3(1.0, 1.0, 10.0));
	return min(da, min(db, dc));

}

// Sponge (slow)
float mSponge(vec3 p) {
	float d = sdBox(p, vec3(1.0));

	float s = 1.0;
	for(int m = 0; m < 4; m++) {
		vec3 a = mod(p * s, 2.0) - 1.0;
		s *= 3.0;
		vec3 r = 1.0 - 3.0 * abs(a);

		float c = sdCross(r) / s;
		d = max(d, c);
	}

	return d;
}

// Sponge
float menger(vec3 pos) {
	vec3 pos2 = abs(pos *= 0.5) - 0.5;
	pos += 0.5;
	float d = max(max(pos2.x, pos2.y), pos2.z);
	float p = 1.0;
	for(int i = 1; i < 4; i++) {
		p *= 3.0;
		pos2 = 0.5 - abs(mod(pos * p, 3.0) - 1.5);
		d = max(min(max(pos2.x, pos2.z), min(max(pos2.x, pos2.y), max(pos2.y, pos2.z))) / p, d);
	}
	return d * 2.0;
}

// Mandelbulb
float mandelbulb2(vec3 p) {
	vec3 w = p;
	float m = dot(w, w);

	float dz = 1.0;

	for(int i = 0; i < 4; i++) {
#if 1
        // polynomial version (no trigonometrics, but MUCH slower)
		float m2 = m * m;
		float m4 = m2 * m2;
		dz = 8.0 * sqrt(m4 * m2 * m) * dz + 1.0;

		float x = w.x;
		float x2 = x * x;
		float x4 = x2 * x2;
		float y = w.y;
		float y2 = y * y;
		float y4 = y2 * y2;
		float z = w.z;
		float z2 = z * z;
		float z4 = z2 * z2;

		float k3 = x2 + z2;
		float k2 = inversesqrt(k3 * k3 * k3 * k3 * k3 * k3 * k3);
		float k1 = x4 + y4 + z4 - 6.0 * y2 * z2 - 6.0 * x2 * y2 + 2.0 * z2 * x2;
		float k4 = x2 - y2 + z2;

		w.x = p.x + 64.0 * x * y * z * (x2 - z2) * k4 * (x4 - 6.0 * x2 * z2 + z4) * k1 * k2;
		w.y = p.y + -16.0 * y2 * k3 * k4 * k4 + k1 * k1;
		w.z = p.z + -8.0 * y * k4 * (x4 * x4 - 28.0 * x4 * x2 * z2 + 70.0 * x4 * z4 - 28.0 * x2 * z2 * z4 + z4 * z4) * k1 * k2;
#else
        // trigonometric version (MUCH faster than polynomial)

        // dz = 8*z^7*dz
		dz = 8.0 * pow(m, 3.5) * dz + 1.0;

        // z = z^8+c
		float r = length(w);
		float b = 8.0 * acos(w.y / r);
		float a = 8.0 * atan(w.x, w.z);
		w = p + pow(r, 8.0) * vec3(sin(b) * sin(a), cos(b), sin(b) * cos(a));
#endif        

        // trap = min( trap, vec4(abs(w),m) );

		m = dot(w, w);
		if(m > 256.0)
			break;
	}

    // distance estimation (through the Hubbard-Douady potential)
	return 0.25 * log(m) * sqrt(m) / dz;
}

float sdSierpinski(vec3 p) {
	int depth = 15;
	vec3 a1 = vec3(0., +2., 0.);
	vec3 a2 = vec3(-1.4, 0., -1.0);
	vec3 a3 = vec3(+1.4, 0., -1.0);
	vec3 a4 = vec3(0., 0., +1.4);
	float scale = 2.;
	for(int iter = 0; iter < depth; ++iter) {
		vec3 c = a1;
		float dist = length(p - a1);
		float d = length(p - a2);
		if(d < dist) {
			c = a2;
			dist = d;
		}
		d = length(p - a3);
		if(d < dist) {
			c = a3;
			dist = d;
		}
		d = length(p - a4);
		if(d < dist) {
			c = a4;
			dist = d;
		}
		p = scale * p - c * (scale - 1.);
	}
	return length(p) / pow(scale, float(depth)) - 0.002;
}

vec2 map(vec3 p) {
	// Plane
	vec2 plane = vec2(sdPlane(p, vec3(0.0, 1.0, 0.0), 1.0), 1.0);

	// Sphere + Box
	vec3 sphere_pos = vec3(2.0, 0.0, 0.0);
	vec2 sphere = vec2(sdSphere(p - sphere_pos, 0.75), 1.0);
	vec2 box = vec2(sdRoundBox(p - sphere_pos + vec3(0.0, 0.7, 0.0), vec3(1.0, 0.4, 1.0), 0.3), 3.0);
	float db = sin(60.0 * p.x) * sin(60.0 * p.y) * sin(60.0 * p.z);
	box.x += db * 0.003;
	vec2 stand = vec2(opSS(sphere.x, box.x, 0.5), 5.0);

	// Sponge
	vec3 sponge_pos = vec3(-5.5, 0.25, 0.5);
	// // vec2 sponge = vec2(mSponge(p - sponge_pos), 4.0f);
	vec2 sponge = vec2(menger((p - sponge_pos) * 0.8) * 1.25, 7.0);

	// Sierpinski
	vec3 sierpinski_pos = vec3(-1.6, -1.0, 0.0);
	vec2 sierpinski = vec2(sdSierpinski((p - sierpinski_pos)), 7.0);

	// Mandelbulb
	vec3 mb_pos = vec3(-5.0, 1.6, 6.0);
	// // vec3 mb_pos = vec3(0.0f, 0.0f, 0.0f);
	vec2 mb = vec2(mandelbulb2((p - mb_pos) * 0.4) * 2.5, 2.0);

	// Glass
	vec3 glass1_pos = vec3(0.0, 0.8, 2.0);
	vec2 glass_ball = vec2(sdSphere(p - glass1_pos, 0.9), 9.0);
	float an = 1.8;
	vec2 c = vec2(sin(an), cos(an));
	vec2 glass_ring = vec2(sdCappedTorus(p - glass1_pos + vec3(0.0, -0.99, 0.0), c, 0.6, 0.15), 4.0);
	vec2 glass1 = vec2(opUS(glass_ball.x, glass_ring.x, 0.4), 9.0);

	// Cone
	vec3 cone_pos = glass1_pos + vec3(0.0, -0.5, 0.0);
	vec2 cone = vec2(sdCone(p - cone_pos, vec2(0.45, 0.25), 1.5), 4.0);


	// Glass 2
	vec3 glass2_pos = sponge_pos + vec3(0.0, 2.5, 0.0);
	vec2 glass2 = vec2(sdTorus(p - glass2_pos, vec2(0.9, 0.3)), 9.0);

	// Gem
	vec3 gem_pos = vec3(2.0, 0.0, -0.5);
	vec2 gem = vec2(sdSolidAngle(p - gem_pos, (vec2(3., 4.) * 0.2), 0.8), 11.0);

	// Hill1
	vec3 hill1_pos = vec3(3.0, -35.0, 37.0);
	vec2 hill1 = vec2(sdSphere(p - hill1_pos, 40.5), 11.0);
	float d1 = sin(0.2 * p.x) * sin(0.2 * p.y) * sin(0.2 * p.z);
	hill1.x += d1 * 4.0;

	// Hill2
	vec3 hill2_pos = vec3(-35.0, -35.0, 10.0);
	vec2 hill2 = vec2(sdSphere(p - hill2_pos, 40.5), 11.0);
	float d2 = sin(0.18 * p.x) * sin(0.18 * p.y) * sin(0.18 * p.z);
	hill2.x += d2 * 4.0;

	// Hill3
	vec3 hill3_pos = vec3(-20.0, -37.0, 30.0);
	vec2 hill3 = vec2(sdSphere(p - hill3_pos, 40.5), 11.0);
	float d3 = sin(0.15 * p.x) * sin(0.15 * p.y) * sin(0.15 * p.z);
	hill3.x += d3 * 3.0;

	// Hill4
	vec3 hill4_pos = vec3(10.0, -35.0, -40.0);
	vec2 hill4 = vec2(sdSphere(p - hill4_pos, 46.5), 11.0);
	float d4 = sin(0.20 * p.x) * sin(0.20 * p.y) * sin(0.20 * p.z);
	hill4.x += d4 * 4.0;

	// Constructive Solid Geometry
	vec3 cs_box_pos = vec3(2.6, 0.5, 3.0);
	vec2 cs_box = vec2(sdBox(p - cs_box_pos, vec3(0.75)), 3.0);
	vec2 cs_sphere = vec2(sdSphere(p - cs_box_pos, 1.0), 7.0);
	// vec2 cs_sphere2 =  vec2(sdSphere(p - cs_box_pos, 0.85f), 6.0f);
	vec2 cs_cross = vec2(sdCross((p - cs_box_pos) * 2.5) / 2.5, 5.0);
	vec2 cs_figure = opIntersection(cs_box, cs_sphere);
	// cs_figure = opSubtraction(cs_sphere2, cs_figure);
	cs_figure = opSubtraction(cs_cross, cs_figure);
	// return cs_figure;

	// Smooth spheres
	vec3 ss_sphere_pos = vec3(0.6, -1.3, -1.6);	
	vec2 ss_sphere1 = vec2(sdSphere(p - ss_sphere_pos, 0.55), 4.0);
	vec2 ss_sphere2 = vec2(sdSphere(p - ss_sphere_pos + vec3(0.0, -0.6, 0.0), 0.2), 4.0);
	vec2 ss_sphere = vec2(opUS(ss_sphere1.x, ss_sphere2.x, 0.1), 1.0);

	// lollipop
	vec3 lollipop_pos = ss_sphere_pos + vec3(0.0, 0.7, 0.0);
	vec2 lollipop_stem = vec2(sdVerticalCapsule(p - lollipop_pos, 0.8, 0.05), 4.0);
	vec2 lollipop_head = vec2(sdEllipsoid(p - lollipop_pos + vec3(0.0, -0.7, 0.0), vec3(0.3, 0.3, 0.1)),  2.);
	vec2 lollipop = opUnion(lollipop_stem, lollipop_head);

	// Tower
	vec3 tower_pos = vec3(-3.2, -0.5, -2.65);
	// vec2 tower_base = vec2(sdCappedCylinder(p - tower_pos + vec3(0.0f, 0.32f, 0.0f) , 0.25f, 1.0f), 1.0f);
	vec2 tower_base = vec2(sdRoundedCylinder(p - tower_pos + vec3(0.0, 0.6, 0.0), 0.55, 0.2, 0.25), 9.0);
	vec2 stem = vec2(sdRoundedCylinder(p - tower_pos + vec3(0.0, -0.5, 0.0), 0.025, 0.02, 0.7), 9.0);
	vec2 torus0 = vec2(sdTorus(p - tower_pos, vec2(0.5, 0.15)), 11.0);
	vec2 torus1 = vec2(sdTorus(p - tower_pos + vec3(0.0, -0.37, -0.0), vec2(0.35, 0.1)), 3.0);
	vec2 torus2 = vec2(sdTorus(p - tower_pos + vec3(0.0, -0.65, 0.0), vec2(0.25, 0.08)), 5.0);
	vec2 torus3 = vec2(sdTorus(p - tower_pos + vec3(0.0, -0.88, 0.0), vec2(0.15, 0.07)), 6.0);
	vec2 torus4 = vec2(sdTorus(p - tower_pos + vec3(0.0, -1.05, 0.0), vec2(0.1, 0.05)), 7.0);

	// Union
	vec2 res;
	// res = plane;
	res = vec2(opUS(plane.x, hill1.x, 0.9), 1.0);
	res = vec2(opUS(res.x, hill2.x, 0.9), 1.0);
	res = vec2(opUS(res.x, hill3.x, 0.9), 1.0);
	res = vec2(opUS(res.x, hill4.x, 0.9), 1.0);

	res = opUnion(res, sponge); //sponge

	res = opUnion(res, sierpinski); // sierpinski

	res = opUnion(res, mb); //mb

	res = opUnion(res, stand); //US

	// Glass
	res = opUnion(res, cone);
	res = opUnion(res, glass1);
	res = opUnion(glass2, res);

	// Tower
	res = opUnion(res, tower_base);
	res = opUnion(res, stem);
	res = opUnion(res, torus0);
	res = opUnion(res, torus1);
	res = opUnion(res, torus2);
	res = opUnion(res, torus3);
	res = opUnion(res, torus4);

	// Gem
	res = opUnion(res, gem);

	// CSG
	res = opUnion(res, cs_figure);

	// smooth sphere
	res = opUnion(res, ss_sphere);

	// lollipop
	res = opUnion(res, lollipop);
	

	return res;
}

vec2 rayMarch(vec3 ro, vec3 rd) { // -> distance, color_ID
	vec2 hit;
	vec2 object = vec2(0.0, 0.0);
	for(int i = 0; i < uMaxSteps; i++) {
		vec3 p = ro + object.x * rd;
		hit = map(p);
		object.x += hit.x;
		object.y = hit.y;
		if(abs(hit.x) < uEps)
			break;
		if(object.x > uMaxDist) {
			// object.y = 0.0f;
			break;
		}
	}
	return object;
}

vec3 getNormal(vec3 p) {
	vec2 e = vec2(uEps * 5.0, 0.0);
	vec3 normal = vec3(map(p).x) - vec3(map(p - e.xyy).x, map(p - e.yxy).x, map(p - e.yyx).x);
	return normalize(normal);
}

float getAmbientOcclusion(vec3 p, vec3 normal) {
	float occ = 0.0;
	float weight = 1.0;
	for(int i = 0; i < 8; i++) {
		float len = 0.01 + 0.02 * float(i * i);
		float dist = map(p + normal * len).x;
		occ += (len - dist) * weight;
		weight *= 0.85;
	}
	return 1.0 - clamp(0.65 * occ, 0.0, 1.0);
}

vec3 getLight(vec3 p, vec3 rd, vec3 color) {
	// Light source position
	vec3 lightSource = vec3(-20.0, 40.0, -30.0);

	// Vector from point to light - L
	vec3 L = normalize(lightSource - p);
	
	// Normal of point p - N
	vec3 N = getNormal(p);
	
	// Vector from point to camera - V
	vec3 V = -rd;
	
	// Reflected light - R
	vec3 R = reflect(-L, N);

	vec3 specColor = vec3(0.5);
	vec3 specular = specColor * pow(clamp(dot(R, V), 0.0, 1.0), 10.0);
	vec3 diffuse = color * clamp(dot(L, N), 0.0, 1.0);
	vec3 ambient = color * 0.05;
	vec3 fresnel = 0.10 * color * pow(1.0 + dot(rd, N), 3.0);

	//shadow
	float d = rayMarch(p + N * 0.02, L).x;

	float AO = getAmbientOcclusion(p, N);
	// float AO = 1.f;

	// return fresnel;
	// return vec3(AO*0.9f);
	if(d < length(lightSource - p))
		return (ambient + fresnel) * AO;
	return diffuse + (ambient + specular + fresnel) * AO;

}

vec3 getMaterial(vec3 p, float id) {
	vec3 material_color;
	switch(int(id)) {
	// Plane
		case 1:
			material_color = vec3(0.2 + 0.4 * mod(floor(p.x) + floor(p.z), 2.0));
			break;
	// pink+++
		case 2:

			material_color = vec3(0.32, 0.51, 0.08);
			break;
	//Red
		case 3:
			material_color = vec3(0.8, 0.152941, 0.03921);
			break;

	// Gray
		case 4:
			material_color = vec3(0.6, 0.6, 0.6);
			break;

	// Orange
		case 5:
			material_color = vec3(0.9098, 0.3490, 0.04313);
			break;

	// orange---
		case 6:
			material_color = vec3(0.68, 0.5607, 0.0);
			break;

	// yellow
		case 7:
			material_color = vec3(1.0, 0.7, 0.0);
			break;

	// pink
		case 8:
			material_color = vec3(0.9, 0.4, 0.6);
			break;

	// glass
		case 9:
			material_color = vec3(0.0);
			break;

	// red glass
		case 10:
			material_color = vec3(0.5, 0.0, 0.0);
			break;

		case 11:
			material_color = vec3(1.0, 0.1, 0.0);
			break;
	// Purple
		case 12:
			material_color = vec3(0.5, 0.0, 0.5);
			break;
		default:
			material_color = vec3(0.0, 0.0, 0.0);
			break;
	}
	return material_color;

}

mat3 getCam(vec3 ro, vec3 lookAt) {
	vec3 camF = normalize(vec3(lookAt - ro));
	vec3 camR = normalize(cross(vec3(0, 1, 0), camF));
	vec3 camU = cross(camF, camR);
	return mat3(camR, camU, camF);
}

vec3 render(vec2 uv) {
	vec3 col = vec3(0.0, 0.0, 0.0);

	// Camera position
	vec3 ro = uCamPos;

	vec3 lookAt = uLookAt;
	vec3 rd = getCam(ro, lookAt) * normalize(vec3(uv, uFov));

	vec2 object = rayMarch(ro, rd);

	vec3 Sky = vec3(0.2, 0.2, 0.3);
	if(object.x < uMaxDist) {
		vec3 p = ro + object.x * rd;
		vec3 material = getMaterial(p, object.y);
		col += getLight(p, rd, material);

		// reflection
		if(uReflect == 1) {
			vec3 N = getNormal(p);
			vec3 reflectDir = reflect(rd, N);
			vec2 hit2 = rayMarch(p + reflectDir * uEps * 3., reflectDir);
			if(hit2.x > uEps) {
				if(hit2.x < uMaxDist) {
					vec3 p2 = p + hit2.x * reflectDir;
					vec3 material2 = getMaterial(p2, hit2.y);
					col += getLight(p2, reflectDir, material2) * 0.3;
				} else {
					col += Sky * 0.1;
				}
			}
		}

		// Glass
		vec3 p_new = p;
		vec2 obj_new = object;
		vec3 rd_new = rd;
		const int max_bounce = 3;
		float factor = 1.0;
		for(int i = 0; i < max_bounce; i++) {
			if(obj_new.y >= 9.0 && obj_new.y <= 10.0) {
				vec3 N = getNormal(p_new);
				rd_new = reflect(rd_new, N);
				obj_new = rayMarch(p_new + rd_new * uEps * 3., rd_new);
				if(obj_new.x < uMaxDist) {
					p_new += obj_new.x * rd_new;
					vec3 material_new = getMaterial(p_new, obj_new.y);
					col += getLight(p_new, rd_new, material_new) * 0.5 * factor;
				} else {
					col = mix(col, Sky, 0.5);
				}
				factor *= 0.7f;
			}
		}

	} else {
		col += Sky;
	}

	return col;
}

vec2 getUV(vec2 offset) {
	return ((gl_FragCoord.xy + offset + 0.5) - iResolution.xy * 0.5) / iResolution.y;
}

void main() {

	vec3 color;

	if(uAa != 4) {
		// No AA
		vec2 uv = getUV(vec2(0.0));
		color = render(uv);

	} else if(uAa == 4) {
		// AA x4
		vec4 e = vec4(0.125, -0.125, 0.375, -0.375);
		color = render(getUV(e.xz)) + render(getUV(e.yw)) + render(getUV(e.wx)) + render(getUV(e.zy));
		color *= 0.25;
	}

	// Gamma correction
	color = pow(color, vec3(0.4545));

	fragColor = vec4(color, 1.);
}