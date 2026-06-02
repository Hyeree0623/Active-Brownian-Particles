#ifndef BD_NODE_H
#define BD_NODE_H

typedef struct {
    double x;
    double y;
} Particle;

typedef struct Node {
    int name;
    Particle particle;
    struct Node* next;
} Node;

#endif 
