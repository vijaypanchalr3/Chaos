import pygame as pg
from numpy import cos, linspace, meshgrid, sin, pi, zeros
import os
import matplotlib.pyplot as plt
import sys


w02 = 10/200

def AP1(theta1,theta2,phi1,phi2):
    diff = theta1-theta2
    return (-w02*3*sin(theta1)-w02*sin(theta1-2*theta2)-2*sin(diff)*(phi2*phi2+phi1*phi1*cos(diff)))/(3-cos(2*diff))

def AP2(theta1,theta2,phi1,phi2):
    diff = theta1-theta2
    return (2*sin(diff)*(phi1*phi1*2+2*w02*cos(theta1)+phi2*phi2*cos(diff)))/(3-cos(2*diff))




def load_image(image, scale=2):
    fullname = os.path.join("./", image)
    image = pg.image.load(fullname)

    size = image.get_size()
    size = (size[0] * scale, size[1] * scale)
    # image = pg.transform.scale(image, size)
    # image = image.convert()
    
    return image
 

class DoublePendulum:
    """

    """
    def __init__(self,K11,K22,theta1,theta2,color="#000000"):
        self.mass = 1000
        self.length = 200
        self.theta1 = theta1
        self.theta2 = theta2
        self.phi1 = 0
        self.phi2 = 0
        self.K11 = K11
        self.K22 = K22
        self.h = 0.080
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
    """

    """
    def __init__(self,AP1,AP2):
        pg.init()
        self.window = pg.display.set_mode((1300,730),pg.RESIZABLE)
        self.pendulum1 = DoublePendulum(AP1,AP2,0,pi+0.1,"#000000")
        self.pendulum2 = DoublePendulum(AP1,AP2,0,pi+0.10001,"#333333")
        self.bg = "#ffffff"
        self.fg = "#000000"
        self.ff = pg.font.Font("Lato-BoldItalic.ttf",32)

    def menu(self):
        run = True
        clock = pg.time.Clock()
        
        
        while run:
            for event in pg.event.get():
                if event.type == pg.QUIT:
                    run = False
                    break

            clock.tick(60)
            self.window.fill(self.bg)
            size= self.window.get_size()
                
            # title
            menu_title = self.ff.render("MENU",True,self.bg,self.fg)
            menu_title_size = menu_title.get_size()
            menu_title_rect = pg.draw.rect(self.window,self.fg,(size[0]//2-menu_title_size[0],10,menu_title_size[0]+200,menu_title_size[1]+10),border_radius=10)
                
            self.window.blit(menu_title,(menu_title_rect.x+100,menu_title_rect.y+5))

            # option:1
            option1_title = self.ff.render("length1: ",True,self.bg,self.fg)
            option1_title_rect = pg.draw.rect(self.window,self.fg,(size[0]//4,10,200,50),border_radius=5)
            self.window.blit(option1_title,option1_title_rect)
            
            pg.display.flip()
                

        pg.quit()
        sys.exit()
                
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
            self.pendulum2.draw(self.window)

            pg.display.flip()
            self.pendulum1.update()
            self.pendulum2.update()
        pg.quit()


class Graphs(DoublePendulum):
    """
    
    """
    def __init__(self,ap1,ap2):
        self.theta1 = 0
        self.theta2 = pi+0.1
        self.phi1 = 0
        self.phi2 = 0
        super().__init__(ap1,ap2,self.theta1,self.theta2)
    def streamlines(self):
        x1,x2 = meshgrid(self.theta1,self.phi1)
        x1_,x2_ = meshgrid(self.theta2,self.theta2)
        v = self.K11(x1,x1_,x1,x2_)
        plt.figure()
        plt.streamplot(x1,x2,x1,v, color='k', linewidth=0.8,density=1.5, minlength=0.01, arrowsize=0.8,arrowstyle="->")
        plt.title("stream plot of approximated equation")
        plt.xlabel("$\theta$")
        plt.ylabel("$\phi$")
        plt.show()
    def phivstheta(self,datapoints):
        self.h = 0.01
        data = zeros([datapoints])
        time = 0
        x = zeros([10000])
        for i in range(datapoints):
            data[i]=self.phi1
            x[i] = time
            time+=self.h
            self.update()
        plt.plot(x,data)
        plt.show()
        
    def energyvstheta(self,datapoints):
        self.h = 0.01
        data = zeros([datapoints])
        time = 0
        x = zeros([datapoints])
        for i in range(datapoints):
            data[i]= (0.5*self.mass*self.length*self.length*(self.phi1*self.phi1+self.phi2*self.phi2+self.phi1*self.phi2*cos(self.theta1-self.theta2)))+(self.mass*10*self.length*(cos(self.theta1)*2+cos(self.theta2)))
            x[i] = time
            time+=self.h
            self.update()
        plt.plot(x,data)
        plt.show()
        
        

        
if __name__=='__main__':
    # graph = Graphs(AP1,AP2)
    # graph.energyvstheta(10000)
    sim = Simulation(AP1,AP2)
    sim.menu()
    
