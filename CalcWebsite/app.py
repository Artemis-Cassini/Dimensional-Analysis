from flask import Flask, render_template, request, jsonify
import dimensional_analysis_calc_2 as dac

app = Flask(__name__)

@app.route("/")
def home():
    # only conversions & calculator
    return render_template("index.html")

@app.route("/about")
def about():
    return render_template("about.html")

@app.route("/how-to")
def how_to():
    return render_template("how_to.html")

@app.route("/about-me")
def about_me():
    return render_template("about_me.html")

@app.route("/contact")
def contact():
    return render_template("contact.html")

# … your /api/analyze route here …

@app.route("/api/analyze", methods=["POST"])
def analyze():
    data = request.get_json()
    query = data.get("query", "")
    # call your existing function
    # Note: dac.dimensional_analysis returns a string with \n separators
    full_output = dac.dimensional_analysis(query)
    # split into lines for the frontend
    steps = full_output.split("\n")
    return jsonify(steps=steps)

if __name__ == "__main__":
    print("▶️  Starting Flask…")
    app.run(debug=True)