import numpy as np
import turtle
import matplotlib.pyplot as plt
import random
def randomWalkb(length):
    steps = []
    x,y = 0,0
    walkx,walky = [x],[y]
    for i in range(length):
        x = random.random()
        y = random.random()
        walkx.append(x)
        walky.append(y)
    return [walkx,walky]

walk = randomWalkb(25)
print walk


turtle.speed('slowest')

walk = randomWalkb(25)

for x, y in zip(*walk):
    #multiply by 10, since 1 pixel differences are hard to see
    turtle.goto(x*100,y*100)

turtle.exitonclick()
