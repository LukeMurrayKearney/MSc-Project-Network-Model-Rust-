use csv::Writer;
use crate::random_graphs::{Output, ResultType};
use serde::{Serialize};
use serde_json;
use std::io::Write;
use std::fs::File;

pub fn outbreak_results_csv(output: Output, result_type: ResultType, path: &str) {
    // serialize the results to json
    let file = std::fs::File::create(path)
        .expect("Failed to create csv");
    let mut writer = Writer::from_writer(file);
    let result = match result_type {
        ResultType::SIR => output.sir,
        ResultType::AvgInfections(_) => output.infections,
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