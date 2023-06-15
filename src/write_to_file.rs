use csv::Writer;
use crate::random_graphs::Network;

pub fn outbreak_results_csv(network: Network) {
    // serialize the results to json
    let file = std::fs::File::create("../../csv/BA_results_10_000.csv")
        .expect("Failed to create csv");
    let mut writer = Writer::from_writer(file);
    for row in network.results.iter() {
        let row_record: Vec<String> = row.into_iter().map(|value| value.to_string()).collect();
        writer.write_record(&row_record).expect("Failed to write to file");
    }
    writer.flush().expect("Failed to flush writer");
}