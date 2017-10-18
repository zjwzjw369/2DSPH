#include <glut.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#define M_PI 3.1415926535898 


typedef struct sim_param_t {
	char* fname;
	int nframes;//number of frames
	int npframe;//steps per frame
	float h;//particle size
	float dt;//time step
	float rho0;//reference density
	float k;//bulk modulus
	float mu;//viscosity
	float g;
}sim_param_t;

int get_params(int argc, char** argv, sim_param_t params);
typedef struct sim_state_t {
	int n;
	float mass;
	float* rho;//densities
	float* x;//position
	float* vh;//velocities (half step)
	float* v;//velocities (full step)
	float* a;//acceleration
}sim_state_t;

sim_param_t params;
sim_state_t* state;
sim_state_t* alloc_state(int n) {
	sim_state_t* s;
	s = (sim_state_t*)malloc(sizeof(sim_state_t));
	s->n = n;
	s->rho = (float*)malloc( n * sizeof(float));
	s->x = (float*)malloc(2 * n * sizeof(float));
	s->vh = (float*)malloc(2 * n * sizeof(float));
	s->v = (float*)malloc(2 * n * sizeof(float));
	s->a = (float*)malloc(2 * n * sizeof(float));
	return s;
}
//void free_state(sim_state_t* s);
void compute_density(sim_state_t* s, sim_param_t* params) {
	int n = s->n;
	float *rho = s->rho;
	const float* x = s->x;
	float h = params->h;
	float h2 = h*h;
	float h8 = (h2*h2)*(h2*h2);
	float C = 4 * s->mass / M_PI / h8;
	memset(rho, 0, n * sizeof(float));
	for (int i = 0; i < n; ++i) {
		rho[i] += 4 * s->mass / M_PI / h2;
		for (int j = i + 1; j < n; ++j) {
			float dx = x[2 * i + 0] - x[2 * j + 0];
			float dy = x[2 * i + 1] - x[2 * j + 1];
			float r2 = dx*dx + dy*dy;
			float z = h2 - r2;
			if (z > 0) {
				float rho_ij = C*z*z*z;
				rho[i] += rho_ij;
				rho[j] += rho_ij;
			}
		}
	}// take advange of the symmetry of the update
}

void compute_accel(sim_state_t* state, sim_param_t* params)
{
	// Unpack basic parameters
	const float h = params->h;
	const float rho0 = params->rho0;
	const float k = params->k;
	const float mu = params->mu;
	const float g = params->g;
	const float mass = state->mass;
	const float h2 = h*h;
	// Unpack system state
	const float*  rho = state->rho;
	const float*  x = state->x;
	const float*  v = state->v;
	float*  a = state->a;
	int n = state->n;
	// Compute density and color
	compute_density(state, params);
	// Start with gravity and surface forces
	for (int i = 0; i < n; ++i) {
		a[2 * i + 0] = 0;
		a[2 * i + 1] = -g;
	}
	// Constants for interaction term
	float C0 = mass / M_PI / ((h2)*(h2));
	float Cp = 15 * k;
	float Cv = -40 * mu;
	// Now compute interaction forces
	for (int i = 0; i < n; ++i) {
		const float rhoi = rho[i];
		for (int j = i + 1; j < n; ++j) {
			float dx = x[2 * i + 0] - x[2 * j + 0];
			float dy = x[2 * i + 1] - x[2 * j + 1];
			float r2 = dx*dx + dy*dy;
			if (r2 < h2) {
				const float rhoj = rho[j];
				float q = sqrt(r2) / h;
				float u = 1 - q;
				float w0 = C0 * u / rhoi / rhoj;
				float wp = w0 * Cp * (rhoi + rhoj - 2 * rho0) * u / q;
				float wv = w0 * Cv;
				float dvx = v[2 * i + 0] - v[2 * j + 0];
				float dvy = v[2 * i + 1] - v[2 * j + 1];
				a[2 * i + 0] += (wp*dx + wv*dvx);
				a[2 * i + 1] += (wp*dy + wv*dvy);
				a[2 * j + 0] -= (wp*dx + wv*dvx);
				a[2 * j + 1] -= (wp*dy + wv*dvy);
			}
		}
	}
}

typedef int(*domain_fun_t)(float, float);
int box_indicator(float x, float y)
{
	return (x < 0.5) && (y < 0.5);
}
int circ_indicator(float x, float y)
{
	float dx = (x - 0.5);
	float dy = (y - 0.3);
	float r2 = dx*dx + dy*dy;
	return (r2 < 0.25*0.25);
}

sim_state_t* place_particles(sim_param_t* param,
	domain_fun_t indicatef)
{
	float h = param->h;
	float hh = h / 1.3;
	// Count mesh points that fall in indicated region.
	int count = 0;
	for (float x = 0; x < 1; x += hh)
		for (float y = 0; y < 1; y += hh)
			count += indicatef(x, y);

	// Populate the particle data structure
	sim_state_t* s = alloc_state(count);
	int p = 0;
	for (float x = 0; x < 1; x += hh) {
		for (float y = 0; y < 1; y += hh) {
			if (indicatef(x, y)) {
				s->x[2 * p + 0] = x;
				s->x[2 * p + 1] = y;
				s->v[2 * p + 0] = 0;
				s->v[2 * p + 1] = 0;
				++p;
			}
		}
	}
	return s;
}

void normalize_mass(sim_state_t* s, sim_param_t* param)
{
	s->mass = 1;
	compute_density(s, param);
	float rho0 = param->rho0;
	float rho2s = 0;
	float rhos = 0;
	for (int i = 0; i < s->n; ++i) {
		rho2s += (s->rho[i])*(s->rho[i]);
		rhos += s->rho[i];
	}
	s->mass *= (rho0*rhos / rho2s);
}
sim_state_t* init_particles(sim_param_t* param)
{
	sim_state_t* s = place_particles(param, circ_indicator);
	normalize_mass(s, param);
	return s;
}

static void damp_reflect(int which, float barrier,
	float* x, float* v, float* vh)
{
	// Coefficient of resitiution
	const float DAMP = 0.75;
	// Ignore degenerate cases
	if (v[which] == 0)
		return;
	// Scale back the distance traveled based on time from collision
	float tbounce = (x[which] - barrier) / v[which];
	x[0] -= v[0] * (1 - DAMP)*tbounce;
	x[1] -= v[1] * (1 - DAMP)*tbounce;
	// Reflect the position and velocity
	x[which] = 2 * barrier - x[which];
	v[which] = -v[which];
	vh[which] = -vh[which];
	// Damp the velocities
	v[0] *= DAMP; vh[0] *= DAMP;
	v[1] *= DAMP; vh[1] *= DAMP;
}

static void reflect_bc(sim_state_t* s)
{
	// Boundaries of the computational domain
	const float XMIN = 0.0;
	const float XMAX = 1.0;
	const float YMIN = 0.0;
	const float YMAX = 1.0;
	float*  vh = s->vh;
	float*  v = s->v;
	float*  x = s->x;
	int n = s->n;
	for (int i = 0; i < n; ++i, x += 2, v += 2, vh += 2) {
		if (x[0] < XMIN) damp_reflect(0, XMIN, x, v, vh);
		if (x[0] > XMAX) damp_reflect(0, XMAX, x, v, vh);
		if (x[1] < YMIN) damp_reflect(1, YMIN, x, v, vh);
		if (x[1] > YMAX) damp_reflect(1, YMAX, x, v, vh);
	}
}

void leapfrog_step(sim_state_t* s, double dt) {
	const float* a = s->a;
	float* vh = s->vh;
	float* v = s->v;
	float* x = s->x;
	int n = s->n;
	for (int i = 0; i < 2 * n; ++i) vh[i] += a[i] * dt;
	for (int i = 0; i < 2 * n; ++i)	v[i] = a[i] * dt / 2;
	for (int i = 0; i < 2 * n; ++i) x[i] += vh[i] * dt;
	reflect_bc(s);
}

void leapfrog_start(sim_state_t* s, double dt)
{
	const float*  a = s->a;
	float*  vh = s->vh;
	float*  v = s->v;
	float*  x = s->x;
	int n = s->n;
	for (int i = 0; i < 2 * n; ++i) vh[i] = v[i] + a[i] * dt / 2;
	for (int i = 0; i < 2 * n; ++i) v[i] += a[i] * dt;
	for (int i = 0; i < 2 * n; ++i) x[i] += vh[i] * dt;
	reflect_bc(s);
}




void check_state(sim_state_t* s)
{
	for (int i = 0; i < s->n; ++i) {
		float xi = s->x[2 * i + 0];
		float yi = s->x[2 * i + 1];
		//assert(xi >= 0 || xi <= 1);
		//assert(yi >= 0 || yi <= 1);
	}
}

void display()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(1.0f, 1.0f, 1.0f);
	for (int i = 0; i < params.npframe; ++i) {
		compute_accel(state, &params);
		leapfrog_step(state, params.dt);
		check_state(state);
	}
	glEnable(GL_POINT_SMOOTH);
	glColor3f(0.2, 0.2, 1.0);
	glPointSize(10.0);
	
	glBegin(GL_POINTS); 
	for (int i = 0; i <state->n;++i) {
		glVertex3f(state->x[i*2], state->x[i * 2+1], 0.0f);
	}
	glEnd();
	glFlush();
	glutSwapBuffers();
}

void idle(void)
{

	display();

}

void init()
{
	params.fname = "output.txt";
	params.nframes = 400;
	params.npframe = 60;
	params.h = 0.02;
	params.dt = 0.00005;
	params.rho0 = 1000;
	params.k = 1000.0;
	params.mu = 0.1;
	params.g = 9.8;
	state = init_particles(&params);
	int nframes = params.nframes;
	int npframe = params.npframe;
	float dt = params.dt;
	int n = state->n;
	compute_accel(state, &params);
	leapfrog_start(state, dt);
	check_state(state);
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-0.5f, 2.0f, -1.0f, 1.0f, -1.0f, 1.0);
}


int main(int argc, char** argv)
{

	
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE);	
	glutInitWindowPosition(100, 100);
	glutInitWindowSize(600, 600);
	glutCreateWindow("HelloOpenGL");
	init();
	glutDisplayFunc(display);
	glutIdleFunc(idle);

	glutMainLoop();

	return 0;
}