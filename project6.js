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
vec3 Shade( Material mtl, vec3 position, vec3 normal, vec3 view )
{
	vec3 color = vec3(0,0,0);
	for ( int i=0; i<NUM_LIGHTS; ++i ) {
		// TO-DO: Check for shadows
		Ray shadowRay;
		shadowRay.dir = normalize( lights[i].position - position );
		shadowRay.pos = position + 1e-3 * shadowRay.dir;
		HitInfo hit;
		if ( IntersectRay( hit, shadowRay ) ) {
			// The shadow ray hits something, so the point is in shadow
			continue;
		} else {
			// TO-DO: If not shadowed, perform shading using the Blinn model
			vec3 L = normalize( lights[i].position - position );
			float NdotL = dot( normal, L );
			vec3 diffuse = mtl.k_d * lights[i].intensity * max( NdotL, 0.0 );
			vec3 H = normalize( L + view );
			vec3 specular = mtl.k_s * lights[i].intensity * pow( max( dot( normal, H ), 0.0 ), mtl.n );
			//color += mtl.k_d * lights[i].intensity;	// change this line
			color += diffuse + specular;	
		}
	}
	return color;
}

// Intersects the given ray with all spheres in the scene
// and updates the given HitInfo using the information of the sphere
// that first intersects with the ray.
// Returns true if an intersection is found.
bool IntersectRay( inout HitInfo hit, Ray ray ) {
    hit.t = 1e30; // Initialize with a very large distance
    bool foundHit = false;

    for (int i = 0; i < NUM_SPHERES; ++i) {
		// TO-DO: Test for ray-sphere intersection
        vec3 sphereToRay = ray.pos - spheres[i].center; 
        float a = dot(ray.dir, ray.dir); // Coefficient of t^2
        float b = 2.0 * dot(sphereToRay, ray.dir); // Coefficient of t
        float c = dot(sphereToRay, sphereToRay) - spheres[i].radius * spheres[i].radius; // Constant term

        // Calculate the discriminant to determine intersections
        float discriminant = b * b - 4.0 * a * c; 

        if (discriminant >= 0.0) { // Potential intersection(s)
            // Calculate the two possible solutions for t
            float t0 = (-b - sqrt(discriminant)) / (2.0 * a);
            float t1 = (-b + sqrt(discriminant)) / (2.0 * a);

            // Find the closest valid intersection
            float t = min(t0, t1); 

			// TO-DO: If intersection is found, update the given HitInfo
            if (t > 0.0 && t < hit.t) { // Valid intersection found
                hit.t = t;
                hit.position = ray.pos + ray.dir * t;
                hit.normal = normalize(hit.position - spheres[i].center); 
                hit.mtl = spheres[i].mtl; 
                foundHit = true;
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
			r.dir = normalize( ray.dir - 2.0 * dot( hit.normal, ray.dir ) * hit.normal ); // Reflect the ray
			r.pos = hit.position + 1e-3 * r.dir;
			
			if ( IntersectRay( h, r ) ) {
				// TO-DO: Hit found, so shade the hit point
				clr += Shade( h.mtl, h.position, h.normal, view );
				// TO-DO: Update the loop variables for tracing the next reflection ray
				hit = h;
				ray = r;
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