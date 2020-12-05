#version 330 core
precision highp float;
layout (location = 0) in vec2 position;

out vec2 coord;

uniform mat4 viewMatrix;

void main(void)
{
	coord = position;
	gl_Position = vec4(coord, 0.0, 1.0);
}