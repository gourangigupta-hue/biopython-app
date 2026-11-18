import matplotlib.pyplot as plt


# Chou-Fasman propensity values
# Format: Amino Acid: (Alpha Helix, Beta Sheet)

propensity = {
    "A": (1.42, 0.83),
    "R": (0.98, 0.93),
    "N": (0.67, 0.89),
    "D": (1.01, 0.54),
    "C": (0.70, 1.19),
    "Q": (1.11, 1.10),
    "E": (1.51, 0.37),
    "G": (0.57, 0.75),
    "H": (1.00, 0.87),
    "I": (1.08, 1.60),
    "L": (1.21, 1.30),
    "K": (1.14, 0.74),
    "M": (1.45, 1.05),
    "F": (1.13, 1.38),
    "P": (0.57, 0.55),
    "S": (0.77, 0.75),
    "T": (0.83, 1.19),
    "W": (1.08, 1.37),
    "Y": (0.69, 1.47),
    "V": (1.06, 1.70)
}


# Function to validate protein sequence
def validate_sequence(sequence):

    valid_amino_acids = "ACDEFGHIKLMNPQRSTVWY"

    if len(sequence) == 0:
        return False

    for amino_acid in sequence:

        if amino_acid not in valid_amino_acids:
            return False

    return True




# Function to predict secondary structure
def predict_structure(sequence):

    n = len(sequence)

    # Store predicted regions
    helix_region = [False] * n
    sheet_region = [False] * n

    # STEP 1: Find alpha helix nucleation regions
    for i in range(n - 5):

        window = sequence[i:i+6]

        count = 0

        for amino_acid in window:
            if propensity[amino_acid][0] > 1.0:
                count += 1

        if count >= 4:

            start = i
            end = i + 6

            # Extend alpha helix to the left
            while start > 0:
                values = [
                    propensity[a][0]
                    for a in sequence[start-1:start+3]
                ]

                if sum(values) / 4 >= 1.0:
                    start -= 1
                else:
                    break

            # Extend alpha helix to the right
            while end < n:
                values = [
                    propensity[a][0]
                    for a in sequence[end-3:end+1]
                ]

                if sum(values) / 4 >= 1.0:
                    end += 1
                else:
                    break

            # Check final alpha helix region
            region = sequence[start:end]

            alpha_avg = sum(
                propensity[a][0] for a in region
            ) / len(region)

            beta_avg = sum(
                propensity[a][1] for a in region
            ) / len(region)

            if len(region) >= 6 and alpha_avg > 1.03 and alpha_avg > beta_avg:
                for j in range(start, end):
                    helix_region[j] = True

    # STEP 2: Find beta sheet nucleation regions
    for i in range(n - 4):

        window = sequence[i:i+5]

        count = 0

        for amino_acid in window:
            if propensity[amino_acid][1] > 1.0:
                count += 1

        if count >= 3:

            start = i
            end = i + 5

            # Extend beta sheet to the left
            while start > 0:
                values = [
                    propensity[a][1]
                    for a in sequence[start-1:start+3]
                ]

                if sum(values) / 4 >= 1.0:
                    start -= 1
                else:
                    break

            # Extend beta sheet to the right
            while end < n:
                values = [
                    propensity[a][1]
                    for a in sequence[end-3:end+1]
                ]

                if sum(values) / 4 >= 1.0:
                    end += 1
                else:
                    break

            # Check final beta sheet region
            region = sequence[start:end]

            alpha_avg = sum(
                propensity[a][0] for a in region
            ) / len(region)

            beta_avg = sum(
                propensity[a][1] for a in region
            ) / len(region)

            if len(region) >= 5 and beta_avg > 1.05 and beta_avg > alpha_avg:
                for j in range(start, end):
                    sheet_region[j] = True

    # STEP 3: Generate final prediction
    prediction = ""

    for i in range(n):

        if helix_region[i] and sheet_region[i]:

            if propensity[sequence[i]][0] > propensity[sequence[i]][1]:
                prediction += "H"
            else:
                prediction += "E"

        elif helix_region[i]:
            prediction += "H"

        elif sheet_region[i]:
            prediction += "E"

        else:
            prediction += "C"

    return prediction

        


    # Function to display prediction summary
def show_summary(result):

    total = len(result)

    h_count = result.count("H")
    e_count = result.count("E")
    c_count = result.count("C")

    print("\n----- Prediction Summary -----")

    print("Alpha Helix:", h_count)
    print("Beta Sheet:", e_count)
    print("Coil:", c_count)

    print("Alpha Helix Percentage:", round(h_count/total*100, 2), "%")
    print("Beta Sheet Percentage:", round(e_count/total*100, 2), "%")
    print("Coil Percentage:", round(c_count/total*100, 2), "%")
    
    # Function to display graph
def show_graph(result):

    h = result.count("H")
    e = result.count("E")
    c = result.count("C")

    structures = ["Alpha Helix", "Beta Sheet", "Coil"]
    counts = [h, e, c]

    plt.bar(structures, counts)

    plt.title("Protein Secondary Structure Prediction")
    plt.xlabel("Secondary Structure")
    plt.ylabel("Number of Amino Acids")

    plt.show()


# Main program
sequence = input("Enter protein sequence: ").strip().upper()

if validate_sequence(sequence):

    result = predict_structure(sequence)

    print("Protein Sequence:", sequence)
    print("Predicted Structure:", result)
    
    show_summary(result)
    show_graph(result)


else:
    print("Invalid protein sequence! Please try again.")
    
