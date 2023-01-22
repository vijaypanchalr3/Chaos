import pygame as pg
from numpy import cos, linspace, meshgrid,sin,pi
import os
import matplotlib.pyplot as plt



w02 = 10
def AP1(theta1,theta2,phi1,phi2):
    diff = theta1-theta2
    return (-w02*3*sin(theta1)-w02*sin(theta1-2*theta2)-2*sin(diff)*(phi2*phi2+phi1*phi1*cos(diff)))/(3-cos(2*diff))

def AP2(theta1,theta2,phi1,phi2):
    diff = theta1-theta2
    return (2*sin(diff)*(phi1*phi1*2+2*w02*cos(theta1)+phi2*phi2*cos(diff)))/(3-cos(2*diff))

def load_image(image, scale=1):
    fullname = os.path.join("./", image)
    image = pg.image.load(fullname)

    size = image.get_size()
    size = (size[0] * scale, size[1] * scale)
    # image = pg.transform.scale(image, size)
    # image = image.convert()
    
    return image
 



class Pendulum:

    def __init__(self,K11,K22,theta1,theta2,color="#000000"):
        self.mass = 1000
        self.length = 200
        self.theta1 = theta1
        self.theta2 = theta2
        self.phi1 = 0
        self.phi2 = 0
        self.K11 = K11
        self.K22 = K22
        self.h = 0.0250
        self.origin = (650,200)
        self.image = load_image("bitmap.png")
        self.color = color

    def update(self):
        """

        Runge kutta methods
        """
        k11 = self.h*self.phi1
        k12 = self.h*self.phi2
        l11 = self.h*self.K11(self.theta1,self.theta2,self.phi1,self.phi2)
        l12 = self.h*self.K22(self.theta1,self.theta2,self.phi1,self.phi2)
        k21 = self.h*(self.phi1+(l11*0.5))
        k22 = self.h*(self.phi2+(l12*0.5))
        l21 = self.h*self.K11(self.theta1+(k11*0.5),self.theta2+(k12*0.5),self.phi1+(l11*0.5),self.phi2+(l12*0.5))
        l22 = self.h*self.K22(self.theta1+(k11*0.5),self.theta2+(k12*0.5),self.phi1+(l11*0.5),self.phi2+(l12*0.5))
        k31 = self.h*(self.phi1+(l21*0.5))
        k32 = self.h*(self.phi2+(l22*0.5))
        l31 = self.h*self.K11(self.theta1+(k21*0.5),self.theta2+(k22*0.5),self.phi1+(l21*0.5),self.phi2+(l22*0.5))
        l32 = self.h*self.K22(self.theta1+(k21*0.5),self.theta2+(k22*0.5),self.phi1+(l21*0.5),self.phi2+(l22*0.5))
        k41 = self.h*(self.phi1+l31)
        k42 = self.h*(self.phi2+l32)
        l41 = self.h*self.K11(self.theta1+k31,self.theta2+k32,self.phi1+l31,self.phi2+l32)
        l42 = self.h*self.K22(self.theta1+k31,self.theta2+k32,self.phi1+l31,self.phi2+l32)
        k_1 = (1/6)*(k11+k41+2*(k21+k31))
        k_2 = (1/6)*(k21+k42+2*(k22+k32))
        l_1 = (1/6)*(l11+l41+2*(l21+l31))
        l_2 = (1/6)*(l12+l42+2*(l22+l32))

        self.theta1+=k_1
        self.theta2+=k_2
        self.phi1+=l_1
        self.phi2+=l_2
        

    def draw(self,window):
        x1 = self.origin[0]+self.length*cos((pi*1.5)-self.theta1)
        y1 = self.origin[1]-self.length*sin((pi*1.5)-self.theta1)
        x2 = x1+self.length*cos((pi*1.5)-self.theta2)
        y2 =y1-self.length*sin((pi*1.5)-self.theta2)
        pg.draw.line(window,self.color,start_pos=self.origin,end_pos=(x1+10,y1+10),width =2)
        pg.draw.line(window,self.color,start_pos=(x1+10,y1+10),end_pos=(x2+10,y2+10),width=2)
        window.blit(self.image,(x1,y1))
        window.blit(self.image,(x2,y2))



class Simulation:

    def __init__(self,AP1,AP2):
        pg.init()
        self.window = pg.display.set_mode((1300,730))
        self.pendulum1 = Pendulum(AP1,AP2,0,1,"#000000")
        # self.pendulum2 = Pendulum(AP1,AP2,0,1.00001,"#333333")
        self.bg = "#ffffff"
    def run(self):
        run = True
        clock = pg.time.Clock()

        while run:
            for event in pg.event.get():
                if event.type==pg.QUIT:
                    run=False
                    break

            clock.tick(60)

            self.window.fill(self.bg)

            self.pendulum1.draw(self.window)
            # self.pendulum2.draw(self.window)

            pg.display.flip()
            self.pendulum1.update()
            # self.pendulum2.update()
        pg.quit()


class Graphs:

    def __init__(self):
        self.theta1 = linspace(-2,2,2000)
        self.theta2 = linspace(-2,2,2000)
        self.phi1 = linspace(-6,6,2000)
        self.phi2 = linspace(-6,6,2000)
        
        
    
    def streamlines(self,ap):
        x1,x2 = meshgrid(self.theta1,self.phi1)
        x1_,x2_ = meshgrid(self.theta2,self.theta2)
        v = ap(x1,x1_,x1,x2_)
        plt.figure()
        plt.streamplot(x1,x2,x1,v, color='k', linewidth=0.8,density=1.5, minlength=0.01, arrowsize=0.8,arrowstyle="->")
        plt.title("stream plot of approximated equation")
        plt.xlabel("$\theta$")
        plt.ylabel("$\phi$")
        plt.show()
        
    

        
if __name__=='__main__':
    sim = Simulation(AP1,AP2)
    sim.run()
