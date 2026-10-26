import random
import math

def calculate_pi(num_samples, side):
    inside_circle = 0

    for _ in range(num_samples):
        x = random.uniform(0, side)
        y = random.uniform(0, side)
        distance = x**2 + y**2

        if distance <= (side) ** 2:
            inside_circle += 1

    pi_estimate = (inside_circle / num_samples) * 4
    return pi_estimate

def am(num, itr):
    total = 0
    for _ in range(itr):
        total += calculate_pi(num, 1)
    return total / itr

def gm(num, itr):
    product = 1
    for _ in range(itr):
        product *= calculate_pi(num, 1)
    return product ** (1 / itr)

def hm(num, itr):
    total = 0
    for _ in range(itr):
        pi_estimate = calculate_pi(num, 1)
        if pi_estimate != 0:
            total += 1 / pi_estimate
    return itr / total


if __name__ == "__main__":
    num_samples = 10000000 # You can change this value for more or fewer samples
    side = 1
    print("real value of pi:     ", math.pi)
    print("arithmetic mean:      ", am(num_samples, 10))
    # print("geometric mean:       ", gm(num_samples, 100), "error:", abs(gm(num_samples, 100) - math.pi))
    # print("harmonic mean:        ", hm(num_samples, 100), "error:", abs(hm(num_samples, 100) - math.pi))
    # print("estimated value of pi:", calculate_pi(num_samples, side), "error:", abs(calculate_pi(num_samples, side) - math.pi))