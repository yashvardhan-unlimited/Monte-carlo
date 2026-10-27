# we will bw making sample space of the total we get after we roll two dices.

import random

def roll_two_dice():
    die1 = random.randint(1, 6)
    die2 = random.randint(1, 6)
    return die1 + die2

sample_space = {
    2 : 0,
    3 : 0,
    4 : 0,
    5 : 0,
    6 : 0,
    7 : 0,
    8 : 0,
    9 : 0,
    10 : 0,
    11 : 0,
    12 : 0
}

total_rolls = 360000

real_sample_space = {
    2 : 1,
    3 : 2,
    4 : 3,
    5 : 4,
    6 : 5,
    7 : 6,
    8 : 5,
    9 : 4,
    10 : 3,
    11 : 2,
    12 : 1
}

def create_sample_space():
    for key in sample_space:
        sample_space[key] = 0

    for _ in range(total_rolls):
        rs = roll_two_dice()
        sample_space[rs] += 1

create_sample_space()
for total, count in sample_space.items():
    probability = count / total_rolls
    real_probability = real_sample_space[total] / 36
    print(f"Total: {total}, Probability: {probability:.4f}, Real Probability: {real_probability:.4f}")


def error():
    lst = []
    for k in sample_space:
        lst.append(abs((sample_space[k] / total_rolls) - (real_sample_space[k] / 36)))
    return sum(lst)

print("\n============Error============")
print("Total error in probability when total_rolls are 3600000: ", error())


