import sqlite3

conn = sqlite3.connect('chepro.db')
cursor = conn.cursor()


# Makes a new table for each cell line and add the corresponding data to it
cursor.execute(
    'CREATE TABLE "MCF-7" AS SELECT * FROM Observation WHERE "cell_line" = "MCF-7"')
cursor.execute(
    'CREATE TABLE "NTERA-2 clone D1" AS SELECT * FROM Observation WHERE "cell_line" = "NTERA-2 clone D1"')
cursor.execute(
    'CREATE TABLE "HL-60" AS SELECT * FROM Observation WHERE "cell_line" = "HL-60"')

# drop the cell_line column
cursor.execute('ALTER TABLE MCF-7 DROP COLUMN "cell_line"')
cursor.execute('ALTER TABLE "NTERA-2 clone D1" DROP COLUMN "cell_line"')
cursor.execute('ALTER TABLE HL-60 DROP COLUMN "cell_line"')

conn.commit()
cursor.close()
conn.close()
