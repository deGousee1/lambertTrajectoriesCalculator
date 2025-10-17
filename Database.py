import mysql.connector
from mysql.connector import Error


def fetch_data_from_mysql():
    """Connects to MySQL and fetches all data from a table."""
    global cursor
    connection = None
    try:
        # Establish the connection to the MySQL server
        connection = mysql.connector.connect(
            host='localhost',  # Your MySQL host
            database='TrajectoriesDatabase',  # Your database name
            user='Finix',  # Your username
            password='Lambert123!@'  # Your password
        )

        if connection.is_connected():
            cursor = connection.cursor()

            # Execute a SELECT query
            query = "SELECT name FROM celestial_bodies ORDER BY id;"
            cursor.execute(query)
            records = cursor.fetchall()
            print("Records fetched successfully:")
            for row in records:
                print(row)

            query = "SELECT name, radius FROM celestial_bodies WHERE name = 'Earth';"
            cursor.execute(query)

            # Fetch all the records from the query result
            records = cursor.fetchall()

            print("Records fetched successfully:")
            for row in records:
                print(row)

    except Error as e:
        print(f"Error while connecting to MySQL: {e}")

    finally:
        # Close the connection
        if connection and connection.is_connected():
            cursor.close()
            connection.close()
            print("MySQL connection is closed.")


if __name__ == '__main__':
    fetch_data_from_mysql()