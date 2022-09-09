# Computer-Graphics-and-Animation
-Example Images: The following images are the outputs of the ray-tracing program using the images from the input files as inputs.
![image](https://raw.githubusercontent.com/shipeng0128/Computer-Graphics-and-Animation/main/images/balls.png) ![image](https://raw.githubusercontent.com/shipeng0128/Computer-Graphics-and-Animation/main/images/gears.png)



-Images with problems: The following images shows the one of the bugs that I spent most time debuging.
    For this image, there are too much noise in it. And there are some parts of the image that did not get correctly generated.
  
  

Description of my implementation:
    For this project, I implemented a ray tracer, which is a computer rendering teahnique that can realistically simulate the lighting of a scene and an object.
To achieve such object, according to the order of my implemenation, at first, the program traverse through the surface of the object to see what a ray will hit 
on the object. To do this, I implemented the intersect() function, which helps the program to check if the ray hits the triangles and spheres that form the 
object together. After accomplishing this step, the output will be an outline of the image. Then, the next step will be the shaing part. For this part, the 
program will check for the shadows first. Then, the program will do the shading. For the code of shading, most of the logic and code are copied from the first 
assignment of this course.(The shading assignment). Finally, the last part I implemented is the reflection and the refraction.

The problems that I encountered:
    To me, the biggest problem I faced is debugging. Most of the time, it takes me a lot of time to debug. When implementing the intersect() function,
I am having a bug which can not generate any output. In this case, I have no idea what the problem with my code is. In addition, the running time of the 
program troubles me a lot as well. After implementing most parts of the program, I start to debug the program as a whole. At this time, it takes several minutes
for my laptop to generate a output. As a result, it makes my debugging process very inefficient. 
    When writin the code with some complex formula, it is very confusing to me that I have to implement a calculation seperately. For example, when writing the shading, the program can only run with code like this.
![image](https://raw.githubusercontent.com/shipeng0128/Computer-Graphics-and-Animation/main/images/code.jpg)
I used to write these three lines of code into one line. However, in this case, the program can not compile.
    Finally, the last big problem I encountered is the color of the ball. The color of it is always blue. I checked my code over and over again and asked TA Yiwen for help. But, I could not find the problem with it. Also, because of the long running time of my program, the efficiency of debugging is pretty slow. As a result, I failed to solve this problem.
    
What I learned from this project:
    During the last 9 weeks of learning, I have a general understanding of computer visualizaton about what it is doing and some basic concepts behind it. For the ray-tracing assignment, I learned the logic and theory of how to trace a ray and render a image from it. Even thought the topic of this research is Computer Graphics and Animation, I also learned a lot of things that are not directly related to the topic. For example, I learned a lot of concepts about vectors. Because my high school was an international high school, I was not taught about vectors. This project is my first experience with vectors. What's more, I also learned many new debugging techniques. In addition, since my program can not generate a image at first, I have to try to print out some variables using std::cout function. In conclusion, what I have learned from this project is far more than just the concepts of computer visulization. The programming skills and experience that I got from this research project are also very important to me.
    
In the future:
    I think my program of ray-tracing has much room for improvement. To begin with, I will try to speed up my program using a bounding volume hierarchy. The long running time troubles throughout the whole project. It made debugging takes a long time as well. Secondlly, I will solve the color issue of the ball image. With a faster running speed, I believe it will be much easier and more effective for me to debug. 
