top_length = 0
bottom_length = 100

with open("all.50.repeats", "r") as in_file:
        lines = in_file.readlines()
        source = ""
        for l in lines:
                if "SOURCE:" in l:
                        source = l.split(".")[1][1:].rstrip()
                elif "Start" not in l and "r" in l:
                        data = l.split()
                        length = int(data[2])
                        if length > top_length:
                                top_length = length
                        if length < bottom_length:
                                bottom_length = length

print(top_length)
print(bottom_length)