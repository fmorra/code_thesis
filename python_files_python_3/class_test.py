class Robot:
    def __init__(self, name, color, weight, height): # Constructor
        self.name = name # Instance variables
        self.color = color
        self.weight = weight
        self.height = height

    def introduce_self(self):
        print("My name is: " + self.name)

    def robot_bmi(self):
        BMI = self.height/self.weight
        print("My robot BMI is: " + str(BMI))
        return BMI

class Person:
    def __init__(self, name, personality, sitting, robot): # Constructor
        self.name = name # Instance variables
        self.personality = personality
        self.sitting = sitting
        self.robot = robot

    def sit_down(self):
        if self.sitting is False:
            self.sitting = True

    def stand_up(self):
        if self.sitting is True:
            self.sitting = False


r1 = Robot("A", "red", 30, 130) # Creates Robot object
r2 = Robot("B", "blue", 45, 150)
r1.introduce_self()
r2.introduce_self()
BMI = r1.robot_bmi()
BMI_2 = BMI + 3
print(BMI_2)

p1 = Person("Alice", "aggressive", False, r2)
p2 = Person("Bob", "meek", False, r1)
p1.robot.introduce_self()


