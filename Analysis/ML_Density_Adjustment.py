from sklearn.datasets import load_iris
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from joblib import dump, load
import pandas as pd

# Load the dataset
file_path = '../Data/density_adjustments.txt'
data = pd.read_csv(file_path, sep=' ', header=None)

# Assign the columns
data.columns = ['Adjusted Density', 'Mean', 'CV', 'Number', 'Set Density', 'Distribution', 'Overlap']

# One-hot encode the 'distribution' column
data_encoded = pd.get_dummies(data, columns=['Distribution'])

# Separate the features and the target variable
X = data_encoded.drop('Set Density', axis=1)  # Drop the target to isolate features
y = data['Set Density']

# Split dataset into training and testing
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=0)

# Initialize and train the Naive Bayes classifier
lr = LinearRegression()
model = lr.fit(X_train, y_train)


# Test the model
predictions = model.predict(X_test)
print(f"Mean Squared Error: {mean_squared_error(y_test, predictions)}")


# Save the model to a file
dump(model, '../Data/naive_bayes_model.joblib')
