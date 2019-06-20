
path = "/home/abetharan/IMPACT/RUNS/kappa"

with open(path + "/fort.12", "r") as file:
    data = file.read().replace('\n', "")


print(data)