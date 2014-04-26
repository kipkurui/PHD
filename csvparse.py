import csv
with open('Chip_Seq-details.csv', 'rb') as csvfile:
    Encode = csv.reader(csvfile, delimiter=',')
    for row in Encode:
        #print row[1]
        if row[3]==" antibody=CTCF":
        #if "CTCF" in row:
            print row
csvfile.close()
