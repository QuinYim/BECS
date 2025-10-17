#include <iostream>
#include <chrono>
#include <cstring>
#include <pbc/pbc.h>
#include <openssl/sha.h>
#include <sstream>
#include <iomanip>

using namespace std;
using namespace std::chrono;

// Define SHA256 digest length
#define SHA256_DIGEST_LENGTH 32

// Define global variables
pbc_param_t param;
pairing_t pairing;
element_t g;            // Generator g
element_t h1, h2;       // Base station public keys h1, h2
element_t b1, b2;       // Base station private keys b1, b2
mpz_t order;            // Group order q

// Initialize system parameters
void system_init() {
    // Initialize parameters using Type A pairing
    pbc_param_init_a_gen(param, 256, 256); // Adjusted for better security

    // Initialize pairing
    pairing_init_pbc_param(pairing, param);

    // Get group order q
    mpz_init(order);
    mpz_set(order, pairing->r);

    // Generate generator g
    element_init_G1(g, pairing);
    element_random(g);

    // Generate base station private keys and public keys
    element_init_Zr(b1, pairing);
    element_init_Zr(b2, pairing);
    element_random(b1);
    element_random(b2);

    element_init_G1(h1, pairing);
    element_init_G1(h2, pairing);
    element_mul_zn(h1, g, b1); // h1 = g^b1
    element_mul_zn(h2, g, b2); // h2 = g^b2
}

// Compute SHA256 hash and map to element
void hash_element(element_t out, const uint8_t *data, int len) {
    uint8_t hash[SHA256_DIGEST_LENGTH];
    SHA256(data, len, hash);
    element_from_hash(out, hash, SHA256_DIGEST_LENGTH);
}

// Device registration
void device_registration(element_t d, element_t &dg, const char *ID_d) {
    // Device selects private key d and computes public key dg
    element_random(d);
    element_init_G1(dg, pairing);
    element_mul_zn(dg, g, d); // dg = g^d

    // Compute H(H(ID_d) || dg)
    element_t hash_ID, hash_total;
    element_init_G1(hash_ID, pairing); // Using G1 for hash representation
    element_init_G1(hash_total, pairing);

    // Hash ID_d
    uint8_t hash_ID_bytes[SHA256_DIGEST_LENGTH];
    SHA256((const unsigned char*)ID_d, strlen(ID_d), hash_ID_bytes);
    hash_element(hash_ID, hash_ID_bytes, SHA256_DIGEST_LENGTH);

    // Get dg bytes
    int dg_len = element_length_in_bytes(dg);
    uint8_t* dg_bytes = new uint8_t[dg_len];
    element_to_bytes(dg_bytes, dg);

    // Concatenate H(ID_d) and dg
    uint8_t* concat = new uint8_t[SHA256_DIGEST_LENGTH + dg_len];
    memcpy(concat, hash_ID_bytes, SHA256_DIGEST_LENGTH);
    memcpy(concat + SHA256_DIGEST_LENGTH, dg_bytes, dg_len);

    // Compute H(H(ID_d) || dg)
    hash_element(hash_total, concat, SHA256_DIGEST_LENGTH + dg_len);

    // Store H(H(ID_d) || dg) (for simplicity, printing it)
    cout << "Device registration complete. Stored hash: ";
    char *hash_str = new char[SHA256_DIGEST_LENGTH * 2 + 1];
    for(int i = 0; i < SHA256_DIGEST_LENGTH; ++i)
        sprintf(hash_str + i*2, "%02x", concat[i]);
    hash_str[SHA256_DIGEST_LENGTH * 2] = '\0';
    cout << hash_str << endl;

    // Clean up
    delete[] dg_bytes;
    delete[] concat;
    delete[] hash_str;
    element_clear(hash_ID);
    element_clear(hash_total);
}

// Pseudonym generation using Schnorr Protocol (Non-Interactive via Fiat-Shamir)
bool pseudonym_generation(element_t d, element_t dg, element_t &x, element_t &y, element_t &Y, element_t &Z) {
    // Step 1: D_d sends (tilde_x = g, tilde_y = dg)
    element_t tilde_x, tilde_y;
    element_init_G1(tilde_x, pairing);
    element_init_G1(tilde_y, pairing);
    element_set(tilde_x, g);
    element_set(tilde_y, dg);

    // Step 2: B_b selects gamma, sends x = gamma * tilde_x
    element_t gamma;
    element_init_Zr(gamma, pairing);
    element_random(gamma);

    element_init_G1(x, pairing);
    element_mul_zn(x, tilde_x, gamma); // x = gamma * tilde_x

    // Step 3: D_d computes y = d * x
    element_init_G1(y, pairing);
    element_mul_zn(y, x, d); // y = d * x

    // Step 4: D_d selects delta, computes Y = delta * x
    element_t delta;
    element_init_Zr(delta, pairing);
    element_random(delta);

    element_init_G1(Y, pairing);
    element_mul_zn(Y, x, delta); // Y = delta * x

    // Step 5: Compute challenge epsilon = H(x || y || Y)
    // Serialize x, y, Y
    int x_len = element_length_in_bytes(x);
    int y_len = element_length_in_bytes(y);
    int Y_len = element_length_in_bytes(Y);
    uint8_t* buffer = new uint8_t[x_len + y_len + Y_len];
    element_to_bytes(buffer, x);
    element_to_bytes(buffer + x_len, y);
    element_to_bytes(buffer + x_len + y_len, Y);

    // Hash to get epsilon
    element_t epsilon;
    element_init_Zr(epsilon, pairing);
    hash_element(epsilon, buffer, x_len + y_len + Y_len);

    delete[] buffer;

    // Step 6: Compute Z = delta + epsilon * d mod q
    element_init_Zr(Z, pairing);
    element_mul(Z, epsilon, d);      // Z = epsilon * d
    element_add(Z, Z, delta);        // Z = delta + Z

    // Step 7: B_b verifies Z * x = Y + epsilon * y
    element_t left, right, temp1, temp2;
    element_init_G1(left, pairing);
    element_init_G1(right, pairing);
    element_init_G1(temp1, pairing);
    element_init_G1(temp2, pairing);

    element_mul_zn(left, x, Z);          // left = Z * x
    element_mul_zn(temp1, y, epsilon);   // temp1 = epsilon * y
    element_add(right, Y, temp1);        // right = Y + temp1

    // Compare left and right
    bool valid = (element_cmp(left, right) == 0);
    if (valid) {
        cout << "Pseudonym generation successful. Pseudonym (x, y) generated." << endl;
    } else {
        cout << "Pseudonym generation failed." << endl;
    }

    // Clean up
    element_clear(tilde_x);
    element_clear(tilde_y);
    element_clear(gamma);
    element_clear(delta);
    element_clear(epsilon);
    element_clear(left);
    element_clear(right);
    element_clear(temp1);
    element_clear(temp2);

    return valid;
}

// Certificate issuance using non-interactive zero-knowledge proof
void certificate_issuance(element_t x, element_t y, element_t &cert_zeta_x, element_t &cert_zeta_y) {
    // Base station uses private keys b1 and b2 to sign x and y
    element_init_G1(cert_zeta_x, pairing);
    element_init_G1(cert_zeta_y, pairing);

    element_mul_zn(cert_zeta_x, x, b1); // cert_zeta_x = b1 * x
    element_mul_zn(cert_zeta_y, y, b2); // cert_zeta_y = b2 * y

    cout << "Certificate issuance complete." << endl;
}

// Identity verification using pairings
bool identity_verification(element_t x, element_t y, element_t cert_zeta_x, element_t cert_zeta_y) {
    // Compute pairings e(cert_zeta_x, g) and e(x, h1)
    element_t pairing1, pairing2;
    element_init_GT(pairing1, pairing);
    element_init_GT(pairing2, pairing);

    pairing_apply(pairing1, cert_zeta_x, g, pairing); // e(cert_zeta_x, g)
    pairing_apply(pairing2, x, h1, pairing);          // e(x, h1)

    // Compute pairings e(cert_zeta_y, g) and e(y, h2)
    element_t pairing3, pairing4;
    element_init_GT(pairing3, pairing);
    element_init_GT(pairing4, pairing);

    pairing_apply(pairing3, cert_zeta_y, g, pairing); // e(cert_zeta_y, g)
    pairing_apply(pairing4, y, h2, pairing);          // e(y, h2)

    // Check if e(cert_zeta_x, g) == e(x, h1) and e(cert_zeta_y, g) == e(y, h2)
    bool verify1 = (element_cmp(pairing1, pairing2) == 0);
    bool verify2 = (element_cmp(pairing3, pairing4) == 0);

    if (verify1 && verify2) {
        cout << "Identity verification successful." << endl;
    } else {
        cout << "Identity verification failed." << endl;
    }

    // Clean up
    element_clear(pairing1);
    element_clear(pairing2);
    element_clear(pairing3);
    element_clear(pairing4);

    return (verify1 && verify2);
}

// Utility function to convert element to hex string
string element_to_hex(element_t e) {
    int len = element_length_in_bytes(e);
    unsigned char* buffer = new unsigned char[len];
    element_to_bytes(buffer, e);
    stringstream ss;
    for(int i = 0; i < len; ++i)
        ss << hex << setw(2) << setfill('0') << (int)buffer[i];
    delete[] buffer;
    return ss.str();
}

// Main function
int main() {
    // Record total start time
    auto total_start = high_resolution_clock::now();

    // System initialization
    auto start = high_resolution_clock::now();
    system_init();
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "System initialization time: " << duration.count() << " microseconds" << endl;

    // Device registration
    element_t d;
    element_init_Zr(d, pairing);
    element_t dg;
    const char *ID_d = "Device123";

    start = high_resolution_clock::now();
    device_registration(d, dg, ID_d);
    end = high_resolution_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout << "Device registration time: " << duration.count() << " microseconds" << endl;

    // Pseudonym generation
    element_t x, y, Y, Z;
    bool pseudonym_valid;

    start = high_resolution_clock::now();
    pseudonym_valid = pseudonym_generation(d, dg, x, y, Y, Z);
    end = high_resolution_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout << "Pseudonym generation time: " << duration.count() << " microseconds" << endl;

    if (!pseudonym_valid) {
        cerr << "Pseudonym generation failed. Exiting." << endl;
        // Clean up
        element_clear(d);
        element_clear(dg);
        element_clear(x);
        element_clear(y);
        element_clear(Y);
        element_clear(Z);
        element_clear(g);
        element_clear(h1);
        element_clear(h2);
        element_clear(b1);
        element_clear(b2);
        mpz_clear(order);
        pairing_clear(pairing);
        pbc_param_clear(param);
        return -1;
    }

    // Certificate issuance
    element_t cert_zeta_x, cert_zeta_y;
    start = high_resolution_clock::now();
    certificate_issuance(x, y, cert_zeta_x, cert_zeta_y);
    end = high_resolution_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout << "Certificate issuance time: " << duration.count() << " microseconds" << endl;

    // Identity verification
    start = high_resolution_clock::now();
    bool verification_result = identity_verification(x, y, cert_zeta_x, cert_zeta_y);
    end = high_resolution_clock::now();
    duration = duration_cast<microseconds>(end - start);
    cout << "Identity verification time: " << duration.count() << " microseconds" << endl;

    if (!verification_result) {
        cerr << "Identity verification failed. Exiting." << endl;
        // Clean up
        element_clear(d);
        element_clear(dg);
        element_clear(x);
        element_clear(y);
        element_clear(Y);
        element_clear(Z);
        element_clear(cert_zeta_x);
        element_clear(cert_zeta_y);
        element_clear(g);
        element_clear(h1);
        element_clear(h2);
        element_clear(b1);
        element_clear(b2);
        mpz_clear(order);
        pairing_clear(pairing);
        pbc_param_clear(param);
        return -1;
    }

    // Record total end time
    auto total_end = high_resolution_clock::now();
    auto total_duration = duration_cast<microseconds>(total_end - total_start);
    cout << "Total execution time: " << total_duration.count() << " microseconds" << endl;

    // Clean up
    element_clear(d);
    element_clear(dg);
    element_clear(x);
    element_clear(y);
    element_clear(Y);
    element_clear(Z);
    element_clear(cert_zeta_x);
    element_clear(cert_zeta_y);
    element_clear(g);
    element_clear(h1);
    element_clear(h2);
    element_clear(b1);
    element_clear(b2);
    mpz_clear(order);
    pairing_clear(pairing);
    pbc_param_clear(param);

    return 0;
}
