struct point {
	GLfloat x;
	GLfloat y;
};

point graph[2000];

for(int i = 0; i < 2000; i++) {
	float x = (i - 1000.0) / 100.0;
	graph[i].x = x;
	graph[i].y = sin(x * 10.0) / (1.0 + x * x);
}
