counter = 0

with open("all.50.repeats", "r") as in_file:
        lines = in_file.readlines()
        source = ""
        for l in lines:
                if "SOURCE:" in l:
                        source = l.split(".")[1][1:].rstrip()
                elif "Start" not in l and "r" in l:
                        data = l.split()
                        data[1] = data[1].strip("r")
                        print(f"{source}\t{data[0]}\t{int(data[0])+int(data[2])}\tir_Pos_{counter}\t0\t+")
                        print(f"{source}\t{int(data[1])-int(data[2])}\t{data[1]}\tir_Neg_{counter}\t0\t-")
                        counter += 1