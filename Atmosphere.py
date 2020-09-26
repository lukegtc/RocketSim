from math import *
class AtmosLayer:
    def __init__(self,alt,a,base_alt,base_temp,base_press,g0,r):
        self.alt = alt
        self.a = a
        self.base_alt = base_alt
        self.base_temp = base_temp
        self.base_press = base_press
        self.g0 = g0
        self.r = r
    def temp(self):
        return self.base_temp + self.a*(self.alt-self.base_alt)
    def press(self):
        if self.a == 0:
            return self.base_press * e**((-self.g0/self.r/self.base_temp)*(self.alt-self.base_alt)*1000)
        else:
            return self.base_press*(((self.base_temp+self.a*(self.alt-self.base_alt))/self.base_temp)**(-self.g0/self.r/(self.a/1000)))


def atmossolver(g0,layers,avals,alt_given,r,base_temp,base_press):


    base_temp1 = [base_temp]
    base_press1 = [base_press]


    i = 0
    while True:


        if alt_given<=0:
            return base_press,base_temp,base_press/r/base_temp

        elif layers[i]<alt_given/1000 <=layers[i+1]:

            from_alt = AtmosLayer(alt_given/1000,avals[i],layers[i],base_temp1[-1],base_press1[-1],g0,r)
            press1 = from_alt.press()
            temp1 = from_alt.temp()
            density1 = press1/from_alt.r/temp1
            break
        elif alt_given/1000>layers[-1]:
            temp2 = 2
            press2 = 2
            return press2,temp2,press2/r/temp2
        else:
            new_alt = AtmosLayer(layers[i+1],avals[i],layers[i],base_temp1[-1],base_press1[-1],g0,r)
            base_press1.append(new_alt.press())
            base_temp1.append(new_alt.temp())
            i+=1
    return press1,temp1,density1