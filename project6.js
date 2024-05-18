var raytraceFS = `
struct Ray {
	vec3 pos;
	vec3 dir;
};

struct Material {
	vec3  k_d;	// diffuse coefficient
	vec3  k_s;	// specular coefficient
	float n;	// specular exponent
};

struct Sphere {
	vec3     center;
	float    radius;
	Material mtl;
};

struct Light {
	vec3 position;
	vec3 intensity;
};

struct HitInfo {
	float    t;
	vec3     position;
	vec3     normal;
	Material mtl;
};

uniform Sphere spheres[ NUM_SPHERES ];
uniform Light  lights [ NUM_LIGHTS  ];
uniform samplerCube envMap;
uniform int bounceLimit;

bool IntersectRay( inout HitInfo hit, Ray ray );

// Shades the given point and returns the computed color.
vec3 Shade(Material mtl, vec3 position, vec3 normal, vec3 view){
	vec3 color = vec3(0,0,0);
	for(int i = 0 ; i < NUM_LIGHTS ; i++){									
		Ray shadowRay;
		shadowRay.pos = position + 1.5*normal;
		shadowRay.dir = normalize(lights[i].position - position);

		HitInfo hit;
		if(!IntersectRay(hit, shadowRay)){
			vec3 lightDir = normalize(lights[i].position - position);
			vec3 h = normalize(lightDir + view);
			float costheta = max(0.0,dot(normal, lightDir));
			float cosfi = max(0.0,dot(normal, h));
			if ((costheta > -1.0) ){
				color += lights[i].intensity*(costheta*mtl.k_d);
				if (cosfi >= 0.0){
					color += lights[i].intensity* mtl.k_s* (cosfi * pow(cosfi,mtl.n));
				}
			}
		}
	}
	return color;
}


//bool Intersect_Sphere(Ray ray, Sphere sphere){}

// Intersects the given ray with all spheres in the scene
// and updates the given HitInfo using the information of the sphere
// that first intersects with the ray.
// Returns true if an intersection is found.
bool IntersectRay( inout HitInfo hit, Ray ray )
{
	hit.t = 1e30;
	bool foundHit = false;
	//float bias = 1e-6;
	//vec3 biased_pos = ray.pos + bias * ray.dir;
	for ( int i=0; i<NUM_SPHERES; ++i ) {
		// TO-DO: Test for ray-sphere intersection
		// TO-DO: If intersection is found, update the given HitInfo
		// Test for ray-sphere intersection 
		bool replaced_t = false;
		ray.dir = normalize(ray.dir);
		vec3 sphereCenter = spheres[i].center;
		float sphereRadius = spheres[i].radius;
		vec3 p_minus_c = ray.pos - sphereCenter;
		float a = dot(ray.dir, ray.dir);
		float b = 2.0 * dot(ray.dir, p_minus_c);
		float c = dot(p_minus_c, p_minus_c) - sphereRadius * sphereRadius;
		float Delta = b*b - 4.0*a*c;
		if (Delta>=0.0){
			float t1 = (-b + sqrt(Delta))/2.0*a;
			float t2 = (-b - sqrt(Delta))/2.0*a;
			if ((t1>=0.0) && (t1<hit.t)){
				hit.t = t1;
				replaced_t = true;
			}
			if ((t2>=0.0) && (t2<hit.t)){
				hit.t = t2;
				replaced_t = true;
			}
			

			if(replaced_t){
				foundHit = true;
				//hit.t = tempt;
				hit.position = ray.pos + hit.t * normalize(ray.dir);
				vec3 new_pos_minus_c = hit.position - sphereCenter;
				hit.normal = normalize(new_pos_minus_c); // normal = ray from center of sphere to position, normalized
				hit.mtl = spheres[i].mtl;
			}
		}
	}
	return foundHit;
}

// Given a ray, returns the shaded color where the ray intersects a sphere.
// If the ray does not hit a sphere, returns the environment color.
vec4 RayTracer( Ray ray )
{
	HitInfo hit;
	if ( IntersectRay( hit, ray ) ) {
		vec3 view = normalize( -ray.dir );
		vec3 clr = Shade( hit.mtl, hit.position, hit.normal, view );
		
		// Compute reflections
		vec3 k_s = hit.mtl.k_s;
		for ( int bounce=0; bounce<MAX_BOUNCES; ++bounce ) {
			if ( bounce >= bounceLimit ) break;
			if ( hit.mtl.k_s.r + hit.mtl.k_s.g + hit.mtl.k_s.b <= 0.0 ) break;
			
			Ray r;	// this is the reflection ray
			HitInfo h;	// reflection hit info
			
			// TO-DO: Initialize the reflection ray
			
			r.pos = hit.position + 1e-2 * hit.normal;
			r.dir = normalize(2.0*dot(view,hit.normal)*hit.normal - view);
			

			if ( IntersectRay( h, r ) ) {
				// TO-DO: Hit found, so shade the hit point
				// TO-DO: Update the loop variables for tracing the next reflection ray
				view = -r.dir;
			
				hit = h;
				
				clr = k_s * Shade(h.mtl, h.position, h.normal, view);
				k_s *= hit.mtl.k_s;
				
			} else {
				// The refleciton ray did not intersect with anything,
				// so we are using the environment color
				clr += k_s * textureCube( envMap, r.dir.xzy ).rgb;
				break;	// no more reflections
			}
		}
		return vec4( clr, 1 );	// return the accumulated color, including the reflections
	} else {
		return vec4( textureCube( envMap, ray.dir.xzy ).rgb, 0 );	// return the environment color
	}
}
`;