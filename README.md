# Computer-Graphics-and-Animation
-Example Images: The following images are the outputs of the ray-tracing program using the images from the input files as inputs.

![image](https://raw.githubusercontent.com/shipeng0128/Computer-Graphics-and-Animation/main/images/balls.png) 
![image](https://raw.githubusercontent.com/shipeng0128/Computer-Graphics-and-Animation/main/images/gears.png)
![image](https://raw.githubusercontent.com/shipeng0128/Computer-Graphics-and-Animation/main/images/tea.png)



-Images with problems: The following images shows the one of the bugs that I spent most time debuging.
    For this image, there are too much noise in it. And there are some parts of the image that did not get correctly generated.
  ![image](https://raw.githubusercontent.com/shipeng0128/Computer-Graphics-and-Animation/main/images/test20.jpg)
  In this image of ball, there many parts of the ball that are not correctly rendered. After checking each part of the code, I found out that is was caused by incorrect implementation of the formula in the intersect() function and the shading part. Besides, the color of this image does not seem right to me. TA Yiwen gave me a lot of suggestions on looking for the problem with it. However, it failed to pinpoint the cause of it and correct it. In the future, I may try to work on this and find the problem.
  

Description of my implementation: 
For this project, I implemented a ray tracer, which is a computer rendering technique that can realistically simulate the lighting of a scene and an object. To achieve such an objective, according to the order of my implementation, at first, the program traverses through the surface of the object to see what a ray will hit on the object. To do this, I implemented the intersect() function, which helps the program check if the ray hits the triangles and spheres that form the object together. After accomplishing this step, the output will be an outline of the image. Then, the next step will be the shading part. For this part, the program will check for the shadows first. Then, the program will do the shading. For the code of shading, most of the logic and code is copied from the first assignment of this course. (The shading assignment). Finally, the last part I implemented is the reflection and the refraction.

The problems that I encountered:
    To me, the biggest problem I faced is debugging. Dubugging takes most time. When implementing the intersect() function, I am having a bug that can not generate any output. In this case, I have no idea what the problem with my code is. In addition, the running time of the program troubles me a lot as well. After implementing most parts of the program, I start to debug the program as a whole. At this time, it takes several minutes for my laptop to generate an output. As a result, it makes my debugging process very inefficient. 
    After finishing implementing all the functions, I started to debug as a whole. It is very complicated for me to debug, because many code may be involved in the bug. I have to rule out all possible causes of the bug one by one. And ray tracer is unlike any other programs that I wrote before. I have no idea how to write some test cases for each part of the code. Therefore, it is very hard for me to debug.
    When writing the code with some complex formula, it is very confusing to me that I have to implement a calculation separately. For example, when writing the shading, the program can only run with code like this.

![image](https://raw.githubusercontent.com/shipeng0128/Computer-Graphics-and-Animation/main/images/code.jpg)

I used to write these three lines of code into one line. However, in this case, the program can not compile.
    Finally, the last big problem I encountered is the color of the ball. The color of it is always blue. I checked my code over and over again and asked TA Yiwen for help. But, I could not find the problem with it. Also, because of the long-running time of my program, the efficiency of debugging is pretty slow. As a result, I failed to solve this problem.
    
What I learned from this project:
    During the last 9 weeks of learning, I have a general understanding of computer visualization what it is doing, and some basic concepts behind it. For the ray-tracing assignment, I learned the logic and theory of how to trace a ray and render an image from it. Even though the topic of this research is Computer Graphics and Animation, I also learned a lot of things that are not directly related to the topic. For example, I learned a lot of concepts about vectors. Because my high school was an international high school, I was not taught about vectors. This project is my first experience with vectors. What's more, I also learned many new debugging techniques. In addition, since my program can not generate an image at first, I have to try to print out some variables using the std::cout function. In conclusion, what I have learned from this project is far more than just the concepts of computer visualization. The programming skills and experience that I got from this research project are also very important to me.
    
In the future:
   I think my program of ray-tracing has much room for improvement. To begin with, I will try to speed up my program using a bounding volume hierarchy. The long-running time troubles me a lot throughout the whole project. It made debugging very time-consuming well. Secondly, I will solve the color issue of the ball image. With a faster running speed, I believe it will be much easier and more effective for me to debug. 
