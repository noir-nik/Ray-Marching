#ifndef SETTINGS_H
#define SETTINGS_H

// ====== Rendering resolution ======
#define WINDOW_WIDTH 960
#define WINDOW_HEIGHT 540

/* ====== Rendering settings for CPU and GPU ====== */
// Anti-aliasing
// == 1 off
// == 4  SSAA x4
const int AA = 4; 

/* ====== Camera settings ====== */
// Camera position
const float POS[3] = {2.0f, 3.0f, -5.0f};

// Camera view direction
const float LOOKAT[3] = {-0.3f, 1.0f, -0.6f};

// Camera FOV
// smaller - wider
const float FOV = 0.8f;


/* Raymarching settings */
const float MAX_DIST = 200.0f;
constexpr int MAX_STEPS = 150;
const float EPS = 0.0001f;

// Reflections for all objects in the scene
const int REFLECT = 0;
#define FRAGMENT_SOURCE "scene.frag"

#endif