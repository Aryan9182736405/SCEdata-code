import pandas as pd

def process_transaction_data(input_file, output_file):
    # Initialize lists to store transaction data
    dates, names, descriptions, withdrawals, additions = [], [], [], [], []

    # Read the input file with utf-8 encoding
    with open(input_file, 'r', encoding='utf-8') as file:
        lines = file.readlines()

    # Process each transaction
    i = 0
    while i < len(lines):
        if lines[i].startswith("You paid") or lines[i].endswith("paid you\n"):
            if lines[i].startswith("You paid"):
                name = lines[i].split("paid")[1].strip()
            else:
                name = lines[i].strip().split("paid you")[0].strip()

            date = lines[i + 1].strip()
            description = lines[i + 2].strip()
            amount = lines[i + 3].strip().replace('$', '').replace(',', '')

            if amount.startswith('-'):
                withdrawal = float(amount[1:])
                addition = 0.0
            elif amount.startswith('+'):
                withdrawal = 0.0
                addition = float(amount[1:])
            else:
                withdrawal = float(amount)
                addition = 0.0

            # Debug statements
            print(f"Name: {name}, Date: {date}, Description: {description}, Withdrawal: {withdrawal}, Addition: {addition}")

            dates.append(date)
            names.append(name)
            descriptions.append(description)
            withdrawals.append(withdrawal)
            additions.append(addition)

            i += 4
        else:
            i += 1

    # Create a DataFrame
    df = pd.DataFrame({
        'Date': dates,
        'Name': names,
        'Description': descriptions,
        'Withdrawal': withdrawals,
        'Addition': additions
    })

    # Save DataFrame to a CSV file
    df.to_csv(output_file, index=False)


input_file = r"C:\Users\aryan\Documents\venmotrans.txt"  # replace with your input file name
output_file = r'C:\Users\aryan\Documents\transactionss.csv'  # replace with your desired output file name
process_transaction_data(input_file, output_file)
