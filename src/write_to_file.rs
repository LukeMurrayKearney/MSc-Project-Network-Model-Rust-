use csv::Writer;
use crate::random_graphs::{Output, ResultType};
use crate::useful_functions::DistributionParameters;
use serde::Serialize;
use serde_json;
use std::io::{Write,Read};
use std::fs::File;

pub fn vector_to_csv(values: Vec<(usize, usize)>, file_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let mut writer = Writer::from_path(file_path)?;
    for row in values {
        writer.write_record(&[row.0.to_string(), row.1.to_string()])?;
    }
    writer.flush()?;
    Ok(())
}

pub fn outbreak_results_csv(output: Output, result_type: ResultType, path: &str) {
    // serialize the results to json
    let file = std::fs::File::create(path)
        .expect("Failed to create csv");
    let mut writer = Writer::from_writer(file);
    let result = match result_type {
        ResultType::SEIR => output.seir,
        ResultType::AvgInfections(_) => output.infections,
        ResultType::SecondaryCases(_) => output.secondary_cases,
    };
    for row in result.iter() {
        let row_record: Vec<String> = row.into_iter().map(|value| value.to_string()).collect();
        writer.write_record(&row_record).expect("Failed to write to file");
    }
    writer.flush().expect("Failed to flush writer");
}

pub fn results_json<T: Serialize>(data: &T, file_path: &str) -> std::io::Result<()> {
    let json_string = serde_json::to_string(data)?;
    let mut file = File::create(file_path)?;
    file.write_all(json_string.as_bytes())?;
    Ok(())
}

fn read_csv_file(file_path: &str) -> Result<Vec<Vec<f64>>, Box<dyn std::error::Error>> {
    let mut file = File::open(file_path)?;
    let mut content = String::new();
    file.read_to_string(&mut content)?;

    let mut reader = csv::ReaderBuilder::new().from_reader(content.as_bytes());
    let mut data: Vec<Vec<f64>> = Vec::new();

    let headers = reader.headers()?.clone();
    data.push(headers.iter().map(|header| header.parse::<f64>().unwrap()).collect());

    for result in reader.records() {
        let record = result?;
        let row: Vec<f64> = record
            .iter()
            .map(|value| value.parse::<f64>())
            .collect::<Result<Vec<f64>, _>>()?;
        data.push(row);
    }

    Ok(data)
}

fn params_json(file_path: &str) -> Result<DistributionParameters, Box<dyn std::error::Error>> {
    let mut file = File::open(file_path)?;
    let mut content = String::new();
    file.read_to_string(&mut content)?;

    let my_struct: DistributionParameters = serde_json::from_str(&content)?;

    Ok(my_struct)
}

pub fn read_rates_mat() -> Vec<Vec<f64>> {
    let file_path = "model_input_files/rates_matrix.csv";
    let rates_mat = match read_csv_file(file_path) {
        Ok(data) => data,
        Err(err) => {
            eprintln!("Error {}", err);
            Vec::new()
        }
    };
    rates_mat
}

pub fn read_params_json(file_path: &str) -> DistributionParameters {
    let file_path = "model_input_files/fitting_parameters1.json";
    let my_struct = match params_json(&file_path) {
        Ok(my_struct) => my_struct,
        Err(err) => {
            eprintln!("Error: {}", err);
            DistributionParameters::new()
        }
    };
    my_struct
}